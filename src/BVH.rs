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
            let mut min = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut max = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for i in 0..3 {
                let point = &triangle[i];
                min.0 = min.0.min(point.x());
                min.1 = min.1.min(point.y());
                min.2 = min.2.min(point.z());
                max.0 = max.0.max(point.x());
                max.1 = max.1.max(point.y());
                max.2 = max.2.max(point.z());
            }
            let mut volume = BoxVolume {
                min: Point::new(min.0, min.1, min.2),
                max: Point::new(max.0, max.1, max.2),
                corresponding: None,
            };
            volume.corresponding = Some(volumes.len());
            volumes.push(volume);
        }
        self.root = Some(Box::new(BVH::build_node(&mut volumes, 0, volumes.len())));
    }

    pub fn build_node(volumes: &mut Vec<BoxVolume>, start: usize, end: usize) -> BVHNode {
        if end - start == 1 {
            return BVHNode::new(volumes[start].clone());
        }
        let mut min = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for i in start..end {
            min.0 = min.0.min(volumes[i].min.x());
            min.1 = min.1.min(volumes[i].min.y());
            min.2 = min.2.min(volumes[i].min.z());
            max.0 = max.0.max(volumes[i].max.x());
            max.1 = max.1.max(volumes[i].max.y());
            max.2 = max.2.max(volumes[i].max.z());
        }
        let axis: usize = []
}