use std::ops::{Index, IndexMut};

use crate::geometric::{Point, Ray, Triangle, Vector};
use crate::light::Color;

#[derive(Clone)]
pub struct BufferItem {
    pub index: usize,
    pub depth: f64,
    pub point: Point,
}

impl Default for BufferItem {
    fn default() -> Self {
        BufferItem {
            index: 0,
            depth: f64::INFINITY,
            point: Point::origin(),
        }
    }
}

pub struct Camera {
    pub picture: Picture<Color>,
    pub buffer: Picture<BufferItem>, // index of triangle, depth
    pub position: Point,
    pub direction: Vector,
    pub up: Vector,
    pub right: Vector,
    pub depth: f64,
    horizontal_half: f64,
    vertical_half: f64,
}

impl Camera {
    pub fn point_at(&mut self, target: &Point) {
        let direction = Vector {
            x: target.x() - self.position.x(),
            y: target.y() - self.position.y(),
            z: target.z() - self.position.z(),
        };
        self.direction = direction.normalize();
        self.right = self.up.cross(self.direction).normalize();
        self.up = self.direction.cross(self.right).normalize();
    }
}

impl Camera {
    pub fn get_ray(&self, p0: u32, p1: u32) -> Ray {
        let x = (p0 as f64 + 0.5) / self.depth;
        let y = (p1 as f64 + 0.5) / self.depth;
        let direction = self.direction
            + self.right * (x - self.horizontal_half)
            + self.up * (self.vertical_half - y);
        Ray {
            start: self.position.clone(),
            direction: direction.normalize(),
        }
    }
    
    pub fn aligned(&self) -> bool {
        self.position == Point::origin() 
            && self.direction == Vector { x: 0.0, y: 0.0, z: 1.0 } 
            && self.up == Vector { x: 0.0, y: 1.0, z: 0.0 }
            && self.right == Vector {x: 1.0, y: 0.0, z: 0.0}
    }
    
    pub fn project(&mut self, triangle: &Triangle, index: usize) {
        assert!(self.aligned());
        // project vertices
        let mut projected_vertices = vec![];
        for i in 0..3 {
            let vertex = &triangle[i];
            let x = vertex.x();
            let y = vertex.y();
            let z = vertex.z();
            let projected = (
                (x / z / self.horizontal_half + 1.0) / 2.0 * self.buffer.width as f64,
                (1.0 - y / z / self.vertical_half) / 2.0 * self.buffer.height as f64,
                z,
            );
            projected_vertices.push(projected);
        }
        // cover buffer
        let min_x = projected_vertices.iter().map(|v| v.0).fold(f64::INFINITY, f64::min) as u32;
        let max_x = projected_vertices.iter().map(|v| v.0).fold(f64::NEG_INFINITY, f64::max) as u32;
        let min_y = projected_vertices.iter().map(|v| v.1).fold(f64::INFINITY, f64::min) as u32;
        let max_y = projected_vertices.iter().map(|v| v.1).fold(f64::NEG_INFINITY, f64::max) as u32;
        for y in min_y.max(0)..=max_y.min(self.picture.height - 1) {
            for x in min_x.max(0)..=max_x.min(self.picture.width - 1) {
                if let Some((depth, point)) = triangle.intersect(&self.get_ray(x, y)) {
                    if depth < self.buffer[(x as usize, y as usize)].depth {
                        self.buffer[(x as usize, y as usize)] = BufferItem { index, depth, point };
                    }
                }
            }
        }
    }
}

pub struct Picture<T> {
    pub width: u32,
    pub height: u32,
    pub pixels: Vec<T>,
}

impl<T> Default for Picture<T> where T: Default + Clone {
    fn default() -> Picture<T> {
        Picture {
            width: 800,
            height: 600,
            pixels: vec![T::default(); 800 * 600],
        }
    }
}

impl<T> Index<(usize, usize)> for Picture<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &T {
        &self.pixels[index.1 * self.width as usize + index.0]
    }
}

impl<T> IndexMut<(usize, usize)> for Picture<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.pixels[index.1 * self.width as usize + index.0]
    }
}

impl Camera {
    pub fn new() -> Camera {
        let picture = Picture::default();
        let depth = picture.width.min(picture.height) as f64;
        Camera {
            horizontal_half: picture.width as f64 / depth / 2.,
            vertical_half: picture.height as f64 / depth / 2.,
            picture,
            buffer: Picture::default(),
            position: Point::new(0.0, 0.0, 0.0),
            direction: Vector { x: 0.0, y: 0.0, z: 1.0 },
            up: Vector { x: 0.0, y: 1.0, z: 0.0 },
            right: Vector { x: 1.0, y: 0.0, z: 0.0 },
            depth,
        }
    }
}