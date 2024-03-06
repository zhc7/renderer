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