extern crate bmp;

use bmp::{Image, Pixel};

use crate::camera::Picture;
use crate::geometric::Point;
use crate::light::Color;
use crate::obj_parser::load;
use crate::object::Properties;
use crate::world::World;

mod geometric;
mod light;
mod object;
mod camera;
mod world;
mod shapes;
mod obj_parser;
mod BVH;


fn save(picture: &Picture<Color>) {
    let mut img = Image::new(picture.width, picture.height);
    for y in 0..picture.height {
        for x in 0..picture.width {
            let color = picture[(x as usize, y as usize)];
            img.set_pixel(x, y, Pixel {
                r: color.r,
                g: color.g,
                b: color.b,
            });
        }
    }
    img.save("output.bmp").unwrap();
}


fn main() {
    // left hind system!!!
    let mut world = World::new();
    // let mut object = shapes::cube();
    let mut object2 = load(r"C:\Users\bjhan\3D Objects\3D Builder\Chess set.obj".as_ref());
    let properties = Properties {
        color: Color::default(),
        ambient: 0.1,
        diffuse: 0.7,
        specular: 0.7,
        shininess: 10.0,
    };
    // object.properties = properties.clone();
    object2.properties = properties.clone();
    // world.add_object(object, Point::new(0.0, 0.0, 10.0));
    world.add_object(object2, Point::new(0.0, -20.0, 150.0));
    world.add_light(light::Light {
        color: light::Color { r: 255, g: 255, b: 255 },
        position: Point::new(5.0, 200.0, 3.0),
    }, Point::new(0.0, 0.0, 0.0));
    world.camera.position = Point::new(-100.0, 60.0, 0.0);
    world.camera.point_at(&Point::new(0.0, 0.0, 150.0));
    world.render();
    save(&world.camera.picture);
}
