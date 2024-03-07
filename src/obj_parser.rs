use std::path::Path;

use crate::geometric::{Point, Triangle};
use crate::object::Object;

pub fn load(path: &Path) -> Object {
    let mut object = Object::new();
    let (models, materials) = tobj::load_obj(
        path,
        &tobj::LoadOptions {
            single_index: true,
            triangulate: true,
            ..Default::default()
        },
    ).unwrap();
    
    let mut base = 0;
    for model in models {
        let mesh = &model.mesh;
        for i in (0..mesh.positions.len()).step_by(3) {
            object.points.push(Point::new(
                mesh.positions[i] as f64,
                mesh.positions[i + 2] as f64,
                mesh.positions[i + 1] as f64,
            ));
        }
        for i in (0..mesh.indices.len()).step_by(3) {
            object.triangles.push(Triangle::new(
                &object.points[mesh.indices[i] as usize + base],
                &object.points[mesh.indices[i + 2] as usize + base],
                &object.points[mesh.indices[i + 1] as usize + base],
            ));
        }
        base = object.points.len();
    }
    
    object
}