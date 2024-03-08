use crate::geometric::{Point, Triangle};
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

pub fn sphere(resolution: usize, radius: f64) -> Object {
    let mut object = Object::new();

    // Generate points
    for i in 1..resolution {
        let theta = (i as f64) * std::f64::consts::PI / (resolution as f64);
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        for j in 0..resolution {
            let phi = (j as f64) * 2.0 * std::f64::consts::PI / (resolution as f64);
            let sin_phi = phi.sin();
            let cos_phi = phi.cos();

            let x = sin_theta * cos_phi * radius;
            let y = sin_theta * sin_phi * radius;
            let z = cos_theta * radius;
            
            let point = Point::new(x, y, z);
            object.points.push(point);
        }
    }
    // Add poles
    object.points.push(Point::new(0.0, 0.0, radius));
    object.points.push(Point::new(0.0, 0.0, -radius));

    // Generate triangles
    for i in 1..resolution - 1 {
        for j in 0..resolution {
            let a = (i - 1) * resolution + j;
            let b = i * resolution + j;
            let c = i * resolution + (j + 1) % resolution;
            let d = (i - 1) * resolution + (j + 1) % resolution;
            object.triangles.push(Triangle::new(
                &object.points[a],
                &object.points[b],
                &object.points[c],
            ));
            object.triangles.push(Triangle::new(
                &object.points[a],
                &object.points[c],
                &object.points[d],
            ));
        }
    }
    // Add top triangles
    for i in 0..resolution {
        object.triangles.push(Triangle::new(
            &object.points[object.points.len() - 2],
            &object.points[i],
            &object.points[(i + 1) % resolution],
        ));
    }
    // Add bottom triangles
    for i in 0..resolution {
        object.triangles.push(Triangle::new(
            &object.points[object.points.len() - 1],
            &object.points[(resolution - 2) * resolution + (i + 1) % resolution],
            &object.points[(resolution - 2) * resolution + i],
        ));
    }
    
    object
}