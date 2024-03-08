use std::collections::BinaryHeap;
use std::fmt::Binary;
use crate::camera::BufferItem;
use crate::geometric::{Point, Ray, Triangle, Vector};

trait Volume {
    fn corresponding(&self) -> Option<usize>;

    fn possibly_intersect(&self, ray: &Ray) -> bool;
    
    fn intersect(&self, ray: &Ray) -> Option<(f64, Point)>;
}

struct SphereVolume {
    center: Point,
    size: f64,
    pub corresponding: Option<usize>,
}

impl SphereVolume {
    fn new(center: Point, size: f64) -> SphereVolume {
        SphereVolume {
            center,
            size,
            corresponding: None,
        }
    }
}

#[derive(Clone)]
struct BoxVolume {
    min: Point,
    max: Point,
    corresponding: Option<usize>,
    triangle: Option<Triangle>,
}

impl BoxVolume {
    fn new(min: Point, max: Point) -> BoxVolume {
        BoxVolume {
            min,
            max,
            corresponding: None,
            triangle: None,
        }
    }
}

impl Volume for BoxVolume {
    fn corresponding(&self) -> Option<usize> {
        self.corresponding
    }

    fn possibly_intersect(&self, ray: &Ray) -> bool {
        let mut t_min = [f64::INFINITY; 3];
        let mut t_max = [f64::NEG_INFINITY; 3];
        for i in 0..3 {
            let inv_d = 1.0 / ray.direction[i];
            let mut t0 = (self.min.index(i) - ray.start.index(i)) * inv_d;
            let mut t1 = (self.max.index(i) - ray.start.index(i)) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            t_min[i] = t_min[i].min(t0);
            t_max[i] = t_max[i].max(t1);
        }
        let t0 = t_min.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
        let t1 = t_max.iter().fold(f64::INFINITY, |a, b| a.min(*b));
        t0 <= t1 && t1 > 0.0
    }

    fn intersect(&self, ray: &Ray) -> Option<(f64, Point)> {
        if let Some(ref triangle) = self.triangle {
            triangle.intersect(ray)
        } else {
            None
        }
    }
}

struct BVHNode {
    volume: Box<dyn Volume>,
    left: Option<Box<BVHNode>>,
    right: Option<Box<BVHNode>>,
}

impl BVHNode {
    fn new(volume: BoxVolume) -> BVHNode {
        BVHNode {
            volume: Box::new(volume),
            left: None,
            right: None,
        }
    }
}

pub struct BVH {
    root: Option<Box<BVHNode>>,
}

impl BVH {
    pub fn new() -> BVH {
        BVH {
            root: None,
        }
    }

    pub fn build(&mut self, triangles: &Vec<&Triangle>) {
        let mut volumes = Vec::new();
        for triangle in triangles {
            let mut min = [f64::INFINITY; 3];
            let mut max = [f64::NEG_INFINITY; 3];
            for i in 0..3 {
                let point = &triangle[i];
                for j in 0..3 {
                    min[j] = min[j].min(point.index(j));
                    max[j] = max[j].max(point.index(j));
                }
            }
            let mut volume = BoxVolume::new(Point::from(min), Point::from(max));
            volume.corresponding = Some(volumes.len());
            volume.triangle = Some((*triangle).clone());
            volumes.push(volume);
        }
        let length = volumes.len();
        self.root = Some(Box::new(BVH::build_node(&mut volumes, 0, length)));
    }

    fn build_node(volumes: &mut Vec<BoxVolume>, start: usize, end: usize) -> BVHNode {
        if end - start == 1 {
            return BVHNode::new(volumes[start].clone());
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for i in start..end {
            let volume = &volumes[i];
            for j in 0..3 {
                min[j] = min[j].min(volume.min.index(j));
                max[j] = max[j].max(volume.max.index(j));
            }
        }
        let axis = Vector::from(max) - Vector::from(min);
        let mut max_axis = 0;
        for i in 1..3 {
            if axis[i] > axis[max_axis] {
                max_axis = i;
            }
        }
        volumes[start..end].sort_by(|a, b| {
            a.min.index(max_axis).partial_cmp(&b.min.index(max_axis)).unwrap()
        });
        let mid = (start + end) / 2;
        let mut node = BVHNode::new(BoxVolume::new(Point::from(min), Point::from(max)));
        node.left = Some(Box::new(BVH::build_node(volumes, start, mid)));
        node.right = Some(Box::new(BVH::build_node(volumes, mid, end)));
        node
    }
    
    pub fn intersect(&self, ray: &Ray) -> BinaryHeap<BufferItem> {
        let mut result = BinaryHeap::new();
        if let Some(ref root) = self.root {
            BVH::intersect_node(root, ray, &mut result);
        }
        result
    }
    
    fn intersect_node(node: &Box<BVHNode>, ray: &Ray, result: &mut BinaryHeap<BufferItem>) {
        if node.volume.possibly_intersect(ray) {
            if let Some(ref left) = node.left {
                BVH::intersect_node(left, ray, result)
            }
            if let Some(ref right) = node.right {
                BVH::intersect_node(right, ray, result)
            }
            if let Some(hit) = node.volume.intersect(ray) {
                if hit.0 > 1e-10 {    // avoid hitting it self
                    if let Some(corresponding) = node.volume.corresponding() {
                        result.push(BufferItem {
                            index: corresponding,
                            depth: hit.0,
                            point: hit.1,
                        });
                    }
                }
            }
        }
    }
}