use crate::world::NormMixMode;


#[derive(Clone, Copy, PartialEq)]
pub enum AntiAliasing {
    SSAA,
    None,
}

#[derive(Clone)]
pub struct RenderArgs {
    pub width: u32,
    pub height: u32,
    pub norm_mode: NormMixMode,
    pub threads: u32,
    pub anti_aliasing: AntiAliasing,
    pub ttl: u32,
    pub camera_depth: Option<f64>,
}

impl Default for RenderArgs {
    fn default() -> RenderArgs {
        RenderArgs {
            width: 800,
            height: 600,
            norm_mode: NormMixMode::Phong,
            threads: 1,
            anti_aliasing: AntiAliasing::None,
            ttl: 5,
            camera_depth: None,
        }
    }
}