import taichi as ti
import os 

ti.init(debug=True, arch=ti.opengl)  # Try to run on opengl

""" allocation """
# constant variable
gravity = ti.Vector.field(2, dtype=float, shape=())
res = 512
quality = 1 # Use a larger value for higher-res simulations
N_particle_material = 3
dt = 1e-5 / quality
N_substeps = int(2e-3 // dt)
brick_width = 0.3
brick_height = 0.2
gravity_coefficient = 150

# particle_material
hardening = 10
jelly_hardening = 0.3
E, nu = 5e3, 0.2  # Young's modulus and Poisson's ratio
mu_0 = E / (2 * (1 + nu)) # Initial Lamé parameters
lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu)) # Initial Lamé parameters

# particles
N_particles = 9000 * quality ** 2
particle_position = ti.Vector.field(2, dtype=ti.f32, shape=N_particles)
particle_velocity = ti.Vector.field(2, dtype=ti.f32, shape=N_particles)
particle_material = ti.field(dtype=int, shape=N_particles)
particle_F = ti.Matrix.field(2, 2, dtype=ti.f32, shape=N_particles) # deformation gradient
particle_C = ti.Matrix.field(2, 2, dtype=ti.f32, shape=(N_particles)) # affine momentum matrix
particle_J = ti.field(dtype=ti.f32, shape=N_particles) # determinant of F, for plastic deformation
pVol = 1.
pRho = 1.
pMass = pVol * pRho


# grid
N_grid = 128 * quality
grid_velocity = ti.Vector.field(2, dtype=ti.f32, shape=(N_grid, N_grid))
grid_mass = ti.field(dtype=ti.f32, shape=(N_grid, N_grid))
dx = 1. / N_grid # grid width, dx
inv_dx = float(N_grid)

# viewer
gui = ti.GUI("MLS-MPM", res=res, background_color=0x198844)
MAX_FRAME = 200000
video_manager = ti.tools.VideoManager(output_dir="./results",
                                framerate=24,
                                automatic_build=False)

# particle_material type 0: fluid 1: jelly 2: snow
def init_scene(mat_type: int = 0):
    @ti.kernel
    def init_particle(mat_type: int): 
        print("mat_type {}".format(mat_type))

        # init particle position
        for i in range(N_particles):
            particle_position[i] = [
                (1. - 2. * ti.random()) * brick_width + 0.5,
                (1. - 2. * ti.random()) * brick_height + 0.5,

            ]

            particle_material[i] = mat_type
            particle_velocity[i] = [0, 0]
            particle_F[i] = ti.Matrix([[1, 0], [0, 1]])
            particle_C[i] = ti.Matrix.zero(float, 2, 2)
            particle_J[i] = 1

    @ti.kernel
    def init_grid():
        for i, j in grid_velocity:
            grid_velocity[i,j] = [0, 0]
            grid_mass[i, j] = 0
        

    init_particle(mat_type)
    init_grid()
    gravity[None] = [0, -1]

def simulate():
    for frame in range(MAX_FRAME):
        handle_interaction()
        render(frame)
        for s in range(N_substeps):
            substep()

def handle_interaction():
    if gui.get_event(ti.GUI.PRESS):
        if gui.event.key == "r":
            print("reset.")
            init_scene(0)
        elif gui.event.key in [ti.GUI.ESCAPE, ti.GUI.EXIT]:
            # save_result()
            exit()
        elif gui.event.key == "l":
            init_scene(0)
        elif gui.event.key == "j":
            init_scene(1)
        elif gui.event.key == "s":
            init_scene(2)
        
def render(frame):
    gui.circles(
        particle_position.to_numpy(),
        radius=1.5,
        palette=[0x9bedf3, 0xfefc38, 0xe5f0e9],
        palette_indices=particle_material,
    )
    img_path = "./results/frames"
    gui.show(os.path.join(img_path, f'{frame:06d}.png'))

def substep():
    reset_grid()
    P2G()
    update_grid()
    G2P()

@ti.kernel
def reset_grid():
    """reset grid velocity and mass every iteration"""
    for i, j in grid_velocity:
        grid_velocity[i, j] = [0, 0]
        grid_mass[i, j] = 0

@ti.kernel
def P2G():
    for p_idx in particle_position:
        cell_idx = (particle_position[p_idx] * inv_dx - 0.5).cast(int)
        # if cell_idx.x < 0 or cell_idx.y < 0:
        #     print(p_idx)
        cell_diff = particle_position[p_idx] * inv_dx - cell_idx.cast(float)
        # Quadratic kernels, [MPM  Eqn. 123]
        weight = [
            0.5 * (1.5 - cell_diff) ** 2,
            0.75 - (cell_diff - 1) ** 2,
            0.5 * (cell_diff - 0.5) ** 2,
        ]
        
        # Compute deformation gradient [MLS_MPM, Eqn 17]
        # F_p^(n+1) =  F_p^(n) + (dt * C_p^[n+1]) * F_p^(n)
        # ? shouldn't update it after G2P
        particle_F[p_idx] = (ti.Matrix.identity(float, 2) + dt * particle_C[p_idx]) @ particle_F[p_idx]

        # Compute current Lamé parameters [MPM Eqn. 86]
        # e = ti.exp(hardening * (1.0 - particle_J[p_idx]))
        e = ti.max(0.1, ti.min(5, ti.exp(10 * (1.0 - particle_J[p_idx]))))


        if particle_material[p_idx] == 1:
            e = jelly_hardening
        mu = mu_0 * e
        lambd = lambda_0 * e
        if particle_material[p_idx] == 0:
            mu = 0.0

        # (svd) Polar decomposition for fixed corotated model
        U, sig, V = ti.svd(particle_F[p_idx])
        # J = ti.math.determinant(particle_F[p_idx])
        J = 1.0
        #?
        for d in ti.static(range(2)):
            new_sig = sig[d, d]
            if particle_material[p_idx] == 2:  # Snow
                new_sig = min(max(sig[d, d], 1 - 2.5e-2), 1 + 4.5e-3)  # Plasticity
            particle_J[p_idx] *= sig[d, d] / new_sig
            sig[d, d] = new_sig
            J *= new_sig
        if particle_material[p_idx] == 0:
            # Reset deformation gradient to avoid numerical instability
            particle_F[p_idx] = ti.Matrix.identity(float, 2) * ti.sqrt(J)
        elif particle_material[p_idx] == 2:
            # Reconstruct elastic deformation gradient after plasticity
            particle_F[p_idx] = U @ sig @ V.transpose()
        
        # first piola stress [MPM Eqn. 52]
        P_stress = 2 * mu * (particle_F[p_idx] - U @ V.transpose()) @ particle_F[p_idx].transpose() + ti.Matrix.identity(float, 2) * lambd * J * (J - 1)

        #? Cauchy stress = (1 / det(F)) * P * F_T [MPM Eqn 38]
        C_stress = (1. / J) * P_stress @ particle_F[p_idx].transpose()

        # force contributed to grid [MLS_MPM Eqn16]
        M_p_inv = 4 * inv_dx * inv_dx
        stress = - dt * pVol * M_p_inv * P_stress
        Q = stress + pMass * particle_C[p_idx]
    

        for i, j in ti.static(ti.ndrange(3, 3)):
            cell_dist = (ti.Vector([i, j]).cast(float) - cell_diff) * dx
            w_ij = weight[i][0] * weight[j][1]

            weightedMass = w_ij * pMass
        
            grid_mass[cell_idx + [i, j]] += weightedMass
            grid_velocity[cell_idx + [i, j]] += weightedMass * particle_velocity[p_idx]

            # fused force contribution and affinie momemtum into Ni(x_p^n)Q_p(x_i - x_p^n) [MLS-MPM equation listed after eqn. 28]
            # where Q_p = dt * M_p^-1 * p.volume * p.cachy_stress + p.mass * p.C
            # instead of our stress, being derived from the energy density, i use the weak form with cauchy stress. converted:
            # p.volume_0 * (dΨ/dF)(Fp)*(Fp_transposed) = p.volume * σ
            # sum up weightedMass * v + weight * Q * dist
            grid_velocity[cell_idx + [i, j]] += w_ij * Q @ cell_dist 


@ti.kernel
def update_grid():
    for i, j in grid_mass:
        if grid_mass[i, j] > 0:
            # convert from momentum to velocity
            grid_velocity[i, j] = grid_velocity[i, j] * (1. / grid_mass[i, j])

            # add gravity 
            grid_velocity[i, j] += dt * gravity[None] * gravity_coefficient

            # handle boundary conditions
            if i < 3 and grid_velocity[i, j][0] < 0:
                grid_velocity[i, j][0] = 0
            if i > (N_grid - 3) and grid_velocity[i, j][0] > 0:
                grid_velocity[i, j][0] = 0
            if j < 3 and grid_velocity[i, j][1] < 0:
                grid_velocity[i, j][1] = 0
            if j > (N_grid - 3) and grid_velocity[i, j][1] > 0:
                grid_velocity[i, j][1] = 0

@ti.kernel
def G2P():
    for p_idx in particle_position:
        cell_idx = (particle_position[p_idx] * inv_dx - 0.5).cast(int)
        cell_diff = particle_position[p_idx] * inv_dx - cell_idx.cast(float)
        weight = [
            0.5 * (1.5 - cell_diff) ** 2,
            0.75 - (cell_diff - 1) ** 2,
            0.5 * (cell_diff - 0.5) ** 2,
        ]
        new_v = ti.Vector.zero(float, 2)
        
        B = ti.Matrix.zero(float, 2, 2)
        for i, j in ti.static(ti.ndrange(3,3)):
            cell_dist = ti.Vector([i, j]).cast(float) - cell_diff
            w_ij = weight[i][0] * weight[j][1]

            # transfer from the grid back to particles, [MPM Eqn 175]
            grid_velocity = grid_velocity[cell_idx + ti.Vector([i, j])]
            new_v += w_ij * grid_velocity

            # constructing affine per-particle momentum matrix C from APIC / MLS-MPM.
            # [APIC page 6] below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
            # where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
            B += w_ij * grid_velocity.outer_product(cell_dist)
        
        # D_inv = 4 * inv_dx * inv_dx
        #?
        D_inv = 4 * inv_dx
        particle_C[p_idx] = B * D_inv
        particle_velocity[p_idx] = new_v 

        # advection
        # particle_position[p_idx] = ti.math.clamp(particle_position[p_idx], [0, 0],[N_grid - 4, N_grid - 4])
        particle_position[p_idx] += dt * particle_velocity[p_idx]
        
        if (particle_position[p_idx].x < 0.001):
            particle_velocity[p_idx].x += 0.001 - particle_position[p_idx].x
        if (particle_position[p_idx].y < 0.001):
            particle_velocity[p_idx].y += 0.001 - particle_position[p_idx].y
        if (particle_position[p_idx].x > 0.999):
            particle_velocity[p_idx].x += 0.999 - particle_position[p_idx].x
        if (particle_position[p_idx].y > 0.999):
            particle_velocity[p_idx].y += 0.999 - particle_position[p_idx].y
            
def save_result():
    pass


if __name__ == "__main__":
    print("[Hint] Press R to reset. Press l/j/s to change particle_material to liquid/jelly/snow")
    init_scene()
    simulate()

