use std::collections::HashMap;
use crate::geometric::{Point, Triangle};
use crate::light::Color;

pub struct Properties {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
}

pub struct Object {
    pub triangles: Vec<Triangle>,
    pub properties: Properties,
}

impl Default for Properties {
    fn default() -> Properties {
        Properties {
            color: Color::default(),
            ambient: 0.0,
            diffuse: 0.0,
            specular: 0.0,
            shininess: 0.0,
        }
    }
}

impl Object {
    pub fn new() -> Object {
        Object {
            triangles: Vec::new(),
            properties: Properties::default(),
        }
    }
    
    pub fn get_points(&mut self) -> Vec<&mut Point> {
        let mut points = Vec::new();
        for triangle in &mut self.triangles {
            for i in 0..3 {
                let point = &mut triangle[i];
                if !points.contains(&point) {
                    points.push(point);
                }
            }
        }
        points
    }

    pub fn calc_triangle_norms(&mut self) {
        for triangle in &mut self.triangles {
            triangle.calc_norm();
        }
    }

    pub fn calc_point_norms(&mut self) {
        let mut weights = HashMap::new();
        let mut normals = HashMap::new();
        for triangle in &self.triangles {
            for i in 0..3 {
                let point = &triangle[i];
                let normal = triangle.normal.unwrap();
                let weight = triangle.size();
                if let Some(n) = normals.get_mut(point) {
                    *n = *n + normal * weight;
                } else {
                    normals.insert(point, normal * weight);
                }
                if let Some(w) = weights.get_mut(point) {
                    *w += weight;
                } else {
                    weights.insert(point, weight);
                }
            }
        }
        for (point, normal) in normals.iter_mut() {
            let weight = weights[point];
            *normal = *normal / weight;
        }
    }
}
