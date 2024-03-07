use crate::geometric::{Point, Ray, Triangle, Vector};


trait Volume {
    fn corresponding(&self) -> Option<usize>;

    fn intersect(&self, ray: &Ray) -> bool;
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

struct BoxVolume {
    min: Point,
    max: Point,
    pub corresponding: Option<usize>,
}

impl BoxVolume {
    fn new(min: Point, max: Point) -> BoxVolume {
        BoxVolume {
            min,
            max,
            corresponding: None,
        }
    }
}

struct BVHNode {
    volume: dyn Volume,
    left: Option<Box<BVHNode>>,
    right: Option<Box<BVHNode>>,
}

impl BVHNode {
    fn new(volume: BoxVolume) -> BVHNode {
        BVHNode {
            volume,
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

    pub fn build(&mut self, triangles: &Vec<Triangle>) {
        let mut volumes = Vec::new();
        for triangle in triangles {
            let mut min = [f64::INFINITY; 3];
            let mut max = [f64::NEG_INFINITY; 3];
            for i in 0..3 {
                let point = &triangle[i];
                for j in 0..3 {
                    min = min.min(point[j]);
                    max = max.max(point[j]);
                }
            }
            let mut volume = BoxVolume::new(Point::from(min), Point::from(max));
            volume.corresponding = Some(volumes.len());
            volumes.push(volume);
        }
        self.root = Some(Box::new(BVH::build_node(&mut volumes, 0, volumes.len())));
    }

    pub fn build_node(volumes: &mut Vec<BoxVolume>, start: usize, end: usize) -> BVHNode {
        if end - start == 1 {
            return BVHNode::new(volumes[start].clone());
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for i in start..end {
            let volume = &volumes[i];
            for j in 0..3 {
                min = min.min(volume.min[j]);
                max = max.max(volume.max[j]);
            }
        }
        
    }
}