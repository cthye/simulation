using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LineRotation : MonoBehaviour
{
    private Vector3 point1, point2, point3;
    public float theta1;
    public float theta2;
    public float l = 5.0f;
    public float m = 5.0f;

    private float w1 = 0.0f; // angular velocity
    private float w2 = 0.0f;
    
    // public float dt = 0.000001f;
    public float dt_inv = 50.0f;


    public float g = 9.8f;
    // private float g = 1f;

    private LineRenderer lr;
    private LineRenderer curve;

    private Vector3[] lrPoints;
    private List<Vector3> curvePoints;

    private bool launched = false;
    
    private void SetLinePoints()
    {
        lrPoints[0] = point1;
        lrPoints[1] = point2;
        lrPoints[2] = point3;

        lr.positionCount = 3;
        lr.SetPositions(lrPoints);
    }

    private void DrawCurve() {
        curvePoints.Add(point3);
        curve.positionCount = curvePoints.Count;
        curve.SetPositions(curvePoints.ToArray());
    }


    void Start()
    {
        theta1 = Mathf.PI / 2.0f;
        theta2 = Mathf.PI / 2.0f;
        point1 = new Vector3(0, 0, 0);
        point2 = point1 + new Vector3(Mathf.Sin(theta1) * l, - Mathf.Cos(theta1) * l, 0);
        point3 = point2 + new Vector3(l * Mathf.Sin(theta2), - l * Mathf.Cos(theta2), 0);

        lr = GetComponent<LineRenderer>();
        lr.SetWidth(0.3f, 0.3f);

        curve = GameObject.Find("curve").GetComponent<LineRenderer>();
        curve.SetWidth(0.01f, 0.1f);

        lrPoints = new Vector3[3];
        curvePoints = new List<Vector3>();

        //Tell the line renderer to update its points
        SetLinePoints();
        // DrawCurve();
    }
    
    void Update()
    {
        if(Input.GetKey("r")) {
            launched = false;
            theta1 = Mathf.PI / 2.0f;
            theta2 = Mathf.PI / 2.0f;
            w1 = 0.0f; // angular velocity
            w2 = 0.0f;
    
            point2 = point1 + new Vector3(Mathf.Sin(theta1) * l, - Mathf.Cos(theta1) * l, 0);
            // point3 = point2 + new Vector3(Mathf.Sin(theta2) * l, - Mathf.Cos(theta1) * l, 0);
            point3 = point2 + new Vector3(l * Mathf.Sin(theta2), - l * Mathf.Cos(theta2), 0);

            SetLinePoints();
            curvePoints.Clear();
            DrawCurve();
            return;
        }

        if(Input.GetKey("l")) {
            launched = true;
        }

        if(launched) {
            // float a_p1 = (- g * 3.0f * m * Mathf.Sin(theta1) - m * g * Mathf.Sin(theta1 - 2.0f * theta2) - 2 * Mathf.Sin(theta1 - theta2) * m * (w2 * w2 * l + w1 * w1 * l * Mathf.Cos(theta1 - theta2)))/ (l * m * (3.0f - Mathf.Cos(2.0f * theta1 - 2.0f * theta2)));
            float a_p1 = - 3.0f / 8.0f * w1 * w2 * Mathf.Sin(theta1 - theta2) - 9.0f * g * Mathf.Sin(theta1) / (8.0f * l) + 3.0f * w2 * Mathf.Sin(theta1 - theta2) / 8.0f;
            // float a_p2 = (2.0f * Mathf.Sin(theta1 - theta2) * ((w1 * w1 * l * 2.0f * m) + g * m * 2.0f * Mathf.Cos(theta1) + w2 * w2 * l * m * Mathf.Cos(theta1 - theta2)))/ (l * m * (3.0f - Mathf.Cos(2.0f * theta1 - 2.0f * theta2)));
            float a_p2 = (3.0f * w1 * w2 * Mathf.Sin(theta1-theta2) - 3.0f * g * Mathf.Sin(theta2) / l + 3 * w1 * Mathf.Sin(theta1 - theta2))/2.0f;

            Debug.Log("p1:" + a_p1);
            Debug.Log("p2:" + a_p2);


            Debug.Log("w1:" + w1);
            Debug.Log("w2:" + w1);
            Debug.Log("det1:" + a_p1 / dt_inv);
            Debug.Log("det2:" + a_p2 / dt_inv);

            w1 += a_p1 / dt_inv;
            w2 += a_p2 / dt_inv;
            Debug.Log("w1:" + w1);
            Debug.Log("w2:" + w2);

            // update theta 
            Debug.Log("theta1:" + theta1);
            Debug.Log("theta2:" + theta2);

            Debug.Log("det1:" + w1 / dt_inv);
            Debug.Log("det2:" + w2 / dt_inv);
            theta1 += w1 / dt_inv;
            theta2 += w2 / dt_inv;
            // if(theta1 > 360.0f) {
            //     float k = theta1 / 360.0f;
            //     theta1 -= 360.0f * k;
            // }
            // if(theta2 > 360.0f) {
            //     float k = theta2 / 360.0f;
            //     theta2 -= 360.0f * k;
            // }

            
            Debug.Log("theta1:" + theta1);
            Debug.Log("theta2:" + theta2);

            point2 = point1 + new Vector3(Mathf.Sin(theta1) * l, - Mathf.Cos(theta1) * l, 0);
            point3 = point2 + new Vector3(l * Mathf.Sin(theta2), - l * Mathf.Cos(theta2), 0);

            SetLinePoints();
            DrawCurve();
        }
    }
}
