use crate::geometric::{Point, Triangle, Vector};
use crate::object::Object;


pub fn cube() -> Object {
    let mut object = Object::new();
    let points = [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ];
    let triangles = [
        [0, 2, 1],
        [0, 3, 2],
        [1, 6, 5],
        [1, 2, 6],
        [5, 7, 4],
        [5, 6, 7],
        [4, 3, 0],
        [4, 7, 3],
        [3, 6, 2],
        [3, 7, 6],
        [4, 1, 5],
        [4, 0, 1],
    ];
    for point in &points {
        object.points.push(Point::new(point[0], point[1], point[2]));
    }
    for triangle in &triangles {
        object.triangles.push(Triangle::new(
            &object.points[triangle[0]],
            &object.points[triangle[1]],
            &object.points[triangle[2]],
        ));
    }
    object
}

pub fn sphere(resolution: usize) -> Object {
    let mut object = Object::new();

    // Generate points
    for i in 0..=resolution {
        let theta = (i as f64) * std::f64::consts::PI / (resolution as f64);
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        for j in 0..=resolution {
            let phi = (j as f64) * 2.0 * std::f64::consts::PI / (resolution as f64);
            let sin_phi = phi.sin();
            let cos_phi = phi.cos();

            let x = sin_theta * cos_phi;
            let y = sin_theta * sin_phi;
            let z = cos_theta;

            object.points.push(Point::new(x, y, z));
        }
    }

    // Generate triangles
    for i in 0..resolution {
        for j in 0..resolution {
            let p1 = i * (resolution + 1) + j;
            let p2 = p1 + (resolution + 1);

            object.triangles.push(Triangle::new(
                &object.points[p1],
                &object.points[p1 + 1],
                &object.points[p2],
            ));
            object.triangles.push(Triangle::new(
                &object.points[p2],
                &object.points[p1 + 1],
                &object.points[p2 + 1],
            ));
        }
    }
    
    // adjust orientation
    for triangle in &mut object.triangles {
        triangle.calc_norm();
        if triangle.normal.unwrap().dot(Vector::from(&triangle.a)) < 0.0 {
            (triangle.a, triangle.c) = (triangle.c.clone(), triangle.a.clone());
        }
    }

    object
}