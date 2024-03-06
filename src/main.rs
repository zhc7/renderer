extern crate bmp;

use bmp::{Image, Pixel};

use crate::camera::Picture;
use crate::geometric::Point;
use crate::world::World;

mod geometric;
mod light;
mod object;
mod camera;
mod world;
mod shapes;


fn save(picture: &Picture) {
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
    let mut object = shapes::cube();
    object.properties.ambient = 0.1;
    object.properties.diffuse = 0.7;
    object.properties.specular = 1.0;
    object.properties.shininess = 1.0;
    world.add_object(object, Point::new(0.0, 0.0, 10.0));
    world.add_light(light::Light {
        color: light::Color { r: 255, g: 255, b: 255 },
        position: Point::new(10.0, 10.0, 3.0),
    }, Point::new(0.0, 0.0, 0.0));
    world.camera.position = Point::new(3.0, 3.0, 0.0);
    world.camera.point_at(&Point::new(0.0, 0.0, 10.0));
    world.render();
    save(&world.camera.picture);
}
