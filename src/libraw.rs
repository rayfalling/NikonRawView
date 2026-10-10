//! LibRaw 的 FFI 绑定与安全封装。
//!
//! # 绑定范围
//!
//! **只绑定访问接口**（`libraw_init` / `libraw_get_*` / `libraw_set_*` /
//! `libraw_dcraw_make_mem_image`），**不复制 `libraw_data_t` 的内部结构体定义**。
//! 那个结构体随 LibRaw 版本变化且体积庞大，在 Rust 侧复刻它意味着升级时会
//! 静默读到错误偏移——这类错误不会崩溃，只会让数值悄悄错掉。
//!
//! 唯一复刻的结构是 [`ProcessedImage`]，它只有 16 字节头且是稳定接口。
//!
//! # 许可
//!
//! LibRaw 为 CDDL-1.0。源码随仓库分发于 `third_party/libraw/`，由 `build.rs` 构建。

use std::ffi::CStr;
use std::os::raw::{c_char, c_int, c_uint, c_void};
use std::path::Path;

/// `libraw_data_t` 在 Rust 侧是不透明指针——我们不依赖它的内部布局。
pub type Handle = *mut c_void;

// ---------------------------------------------------------------------------
// C API
// ---------------------------------------------------------------------------

extern "C" {
    fn libraw_init(flags: c_uint) -> Handle;
    fn libraw_close(lr: Handle);
    /// 窄字符入口。Windows 上改用 `libraw_open_wfile`，故此处仅在非 Windows 使用。
    #[allow(dead_code)]
    fn libraw_open_file(lr: Handle, fname: *const c_char) -> c_int;
    /// Windows 专用：宽字符路径。
    ///
    /// `libraw_open_file` 收的是 `const char*`，在 Windows 上按 ANSI 代码页解释，
    /// **含中文的路径一律打不开**。本项目的照片库路径普遍含中文，因此 Windows 上
    /// 必须走这个入口。
    fn libraw_open_wfile(lr: Handle, fname: *const u16) -> c_int;
    fn libraw_unpack(lr: Handle) -> c_int;
    fn libraw_dcraw_process(lr: Handle) -> c_int;
    fn libraw_dcraw_make_mem_image(lr: Handle, errc: *mut c_int) -> *mut ProcessedImage;
    fn libraw_dcraw_clear_mem(img: *mut ProcessedImage);
    fn libraw_strerror(code: c_int) -> *const c_char;
    /// LibRaw 的 C 符号是驼峰命名，用 `link_name` 映射到 Rust 的蛇形命名
    #[link_name = "libraw_versionNumber"]
    fn libraw_version_number() -> c_int;

    fn libraw_set_output_color(lr: Handle, value: c_int);
    fn libraw_set_output_bps(lr: Handle, value: c_int);
    fn libraw_set_gamma(lr: Handle, index: c_int, value: f32);
    fn libraw_set_no_auto_bright(lr: Handle, value: c_int);
    fn libraw_set_user_mul(lr: Handle, index: c_int, value: f32);
    fn libraw_set_demosaic(lr: Handle, value: c_int);

    fn libraw_get_iwidth(lr: Handle) -> c_int;
    fn libraw_get_iheight(lr: Handle) -> c_int;
    fn libraw_get_raw_width(lr: Handle) -> c_int;
    fn libraw_get_raw_height(lr: Handle) -> c_int;
    fn libraw_get_color_maximum(lr: Handle) -> c_int;
    fn libraw_get_cam_mul(lr: Handle, index: c_int) -> f32;
    fn libraw_get_rgb_cam(lr: Handle, i: c_int, j: c_int) -> f32;
}

/// `libraw_processed_image_t` 的头部。
///
/// 这是 LibRaw 的稳定 C 接口，字段与顺序由上游保证：
/// ```c
/// typedef struct {
///   enum LibRaw_image_formats type;      /* int, 偏移 0  */
///   ushort height, width, colors, bits;  /* 偏移 4,6,8,10 */
///   unsigned int data_size;              /* 偏移 12 */
///   unsigned char data[1];               /* 偏移 16 */
/// } libraw_processed_image_t;
/// ```
#[repr(C)]
pub struct ProcessedImage {
    pub typ: c_int,
    pub height: u16,
    pub width: u16,
    pub colors: u16,
    pub bits: u16,
    pub data_size: u32,
    // 其后的 data[] 不在此结构内声明；用 data_ptr() 取
}

impl ProcessedImage {
    /// 像素数据起点（结构体头部之后）。
    ///
    /// # Safety
    /// 仅当 `self` 指向 `libraw_dcraw_make_mem_image` 的返回值时有效。
    unsafe fn data_ptr(&self) -> *const u8 {
        (self as *const Self as *const u8).add(std::mem::size_of::<Self>())
    }
}

/// 输出色彩空间。取值即 LibRaw / dcraw 的 `-o` 参数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputColor {
    /// 相机原始色彩空间，不做转换。
    Camera = 0,
    Srgb = 1,
    AdobeRgb = 2,
    WideGamut = 3,
    /// ProPhoto（ROMM）。
    ProPhoto = 4,
}

/// 去马赛克算法。取值即 LibRaw / dcraw 的 `-q` 参数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Demosaic {
    Linear = 0,
    Vng = 1,
    Ppg = 2,
    Ahd = 3,
    Dcb = 4,
    /// AHD + 色差中值滤波。
    AhdMedian = 11,
    Dht = 12,
    Aahd = 13,
}

impl Demosaic {
    pub fn name(self) -> &'static str {
        match self {
            Demosaic::Linear => "linear",
            Demosaic::Vng => "vng",
            Demosaic::Ppg => "ppg",
            Demosaic::Ahd => "ahd",
            Demosaic::Dcb => "dcb",
            Demosaic::AhdMedian => "ahd-median",
            Demosaic::Dht => "dht",
            Demosaic::Aahd => "aahd",
        }
    }
}

#[derive(Debug)]
pub enum Error {
    /// 无法读取文件。
    Io(String),
    /// LibRaw 返回了非零错误码。
    LibRaw { code: i32, message: String },
    /// 句柄为空。
    NoHandle,
    /// 输出不是预期的位图格式。
    UnexpectedOutput { typ: i32, colors: u16, bits: u16 },
    /// 路径含内嵌 NUL。
    BadPath(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(s) => write!(f, "无法读取文件：{s}"),
            Error::LibRaw { code, message } => write!(f, "LibRaw 错误 {code}：{message}"),
            Error::NoHandle => write!(f, "libraw_init 返回空句柄"),
            Error::UnexpectedOutput { typ, colors, bits } => {
                write!(f, "输出格式非预期：type={typ} colors={colors} bits={bits}")
            }
            Error::BadPath(s) => write!(f, "路径含内嵌 NUL：{s}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

fn err_message(code: c_int) -> String {
    unsafe {
        let p = libraw_strerror(code);
        if p.is_null() {
            return format!("(无错误描述，code={code})");
        }
        CStr::from_ptr(p).to_string_lossy().into_owned()
    }
}

fn check(code: c_int) -> Result<()> {
    if code == 0 {
        Ok(())
    } else {
        Err(Error::LibRaw { code, message: err_message(code) })
    }
}

/// 打开一个 RAW 文件。
///
/// Windows 上走宽字符入口：`libraw_open_file` 按 ANSI 代码页解释路径，
/// 而本项目的照片库路径含中文（如 `Z:\摄影\Nikon\Z8\...`），用窄字符入口
/// 会一个文件都打不开。
#[cfg(windows)]
fn open_raw(lr: Handle, path: &Path) -> Result<()> {
    use std::os::windows::ffi::OsStrExt;
    let mut wide: Vec<u16> = path.as_os_str().encode_wide().collect();
    wide.push(0);
    unsafe { check(libraw_open_wfile(lr, wide.as_ptr())) }
}

#[cfg(not(windows))]
fn open_raw(lr: Handle, path: &Path) -> Result<()> {
    let c = CString::new(path.to_string_lossy().as_bytes())
        .map_err(|_| Error::BadPath(path.display().to_string()))?;
    unsafe { check(libraw_open_file(lr, c.as_ptr())) }
}

/// LibRaw 的版本号，形如 `0x001600`。
pub fn version_number() -> i32 {
    unsafe { libraw_version_number() }
}

/// 版本号的点分表示。
pub fn version_string() -> String {
    let v = version_number();
    format!("{}.{}.{}", v >> 16, (v >> 8) & 0xFF, v & 0xFF)
}

/// 解码选项。
#[derive(Debug, Clone)]
pub struct Options {
    pub demosaic: Demosaic,
    /// 相机白平衡系数（R/G/B 三个乘数）；`None` 表示用文件中的白平衡。
    pub user_mul: Option<[f32; 3]>,
}

impl Default for Options {
    fn default() -> Self {
        // DHT 在细节与伪色之间较均衡，且对高频细节比 AHD 干净
        Self { demosaic: Demosaic::Dht, user_mul: None }
    }
}

/// 一次解码的结果。
#[derive(Debug, Clone)]
pub struct Decoded {
    /// 输出像素的尺寸。
    ///
    /// **注意不是** `iwidth × iheight`：解码层会按 EXIF 方向旋转输出，竖拍照片
    /// 在这里是 `5520×8280` 而 `iwidth×iheight` 是 `8280×5520`。见 [`Self::rotated`]。
    pub width: usize,
    pub height: usize,
    /// 输出是否相对**未旋转**的有效画幅被转过 90°。
    ///
    /// **裁切边距是针对未旋转画幅定义的**，因此调用方必须据此调整左右与上下，
    /// 否则会切错边。
    pub rotated: bool,
    /// 16 位线性 RGB，行主序，长度 = width × height × 3。
    pub pixels: Vec<u16>,
    /// 相机白平衡系数（R/G/B）。
    pub wb: [f32; 3],
    /// 相机 → 输出空间的 3×3 矩阵（行主序），取自 `libraw_get_rgb_cam`。
    pub rgb_cam: [[f32; 3]; 3],
    pub color_maximum: i32,
    pub demosaic: Demosaic,
    /// 本次解码实际使用的输出色彩空间。
    pub output: OutputColor,
    /// 传感器原始尺寸，含光学黑/掩蔽边框。
    pub raw_width: usize,
    pub raw_height: usize,
}

impl Decoded {
    /// 按索引取像素，越界返回 `None`。
    pub fn at(&self, x: usize, y: usize) -> Option<[u16; 3]> {
        if x >= self.width || y >= self.height {
            return None;
        }
        let i = (y * self.width + x) * 3;
        Some([self.pixels[i], self.pixels[i + 1], self.pixels[i + 2]])
    }

    /// 输出是否比传感器全幅小——即为有效像素区而非含边框的全幅。
    pub fn is_cropped(&self) -> bool {
        self.width < self.raw_width || self.height < self.raw_height
    }
}

/// 一个 LibRaw 句柄的 RAII 包装。
struct Session(Handle);

impl Drop for Session {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { libraw_close(self.0) };
        }
    }
}

/// 把一次 `libraw_dcraw_make_mem_image` 的结果拷成 Rust 侧数据。
///
/// # Safety
/// `img` 必须是当前句柄最近一次 `make_mem_image` 的返回值。
unsafe fn take_image(img: *mut ProcessedImage) -> Result<(usize, usize, Vec<u16>)> {
    let h = &*img;
    if h.colors != 3 || (h.bits != 8 && h.bits != 16) {
        let r = Err(Error::UnexpectedOutput { typ: h.typ, colors: h.colors, bits: h.bits });
        libraw_dcraw_clear_mem(img);
        return r;
    }
    let w = h.width as usize;
    let ht = h.height as usize;
    let n = w * ht * 3;
    let src = h.data_ptr();
    let mut out = vec![0u16; n];
    if h.bits == 16 {
        for (i, o) in out.iter_mut().enumerate() {
            let p = src.add(i * 2);
            *o = u16::from_ne_bytes([*p, *p.add(1)]);
        }
    } else {
        for (i, o) in out.iter_mut().enumerate() {
            *o = (*src.add(i) as u16) * 257; // 8 位按比例扩到 16 位
        }
    }
    libraw_dcraw_clear_mem(img);
    Ok((w, ht, out))
}

/// 解码一张 RAW，以指定输出色彩空间产出 16 位像素，并取出标定元数据。
///
/// 输出恒为线性（gamma 索引 0 与 1 都设为 1.0）且关闭自动亮度，因此结果不含
/// 任何相机内观感或自动增益。
fn decode_with_output(path: &Path, opts: &Options, output: OutputColor) -> Result<Decoded> {
    if !path.is_file() {
        return Err(Error::Io(format!("文件不存在：{}", path.display())));
    }

    let session = Session(unsafe { libraw_init(0) });
    if session.0.is_null() {
        return Err(Error::NoHandle);
    }
    let lr = session.0;

    unsafe {
        open_raw(lr, path)?;
        check(libraw_unpack(lr))?;

        if let Some(m) = opts.user_mul {
            for (i, v) in m.iter().enumerate() {
                libraw_set_user_mul(lr, i as c_int, *v);
            }
        }
        libraw_set_demosaic(lr, opts.demosaic as c_int);

        libraw_set_output_color(lr, output as c_int);
        libraw_set_output_bps(lr, 16);
        libraw_set_gamma(lr, 0, 1.0);
        libraw_set_gamma(lr, 1, 1.0);
        libraw_set_no_auto_bright(lr, 1);
        check(libraw_dcraw_process(lr))?;

        // rgb_cam 在 dcraw_process 内部计算，因此必须在处理后取
        let mut rgb_cam = [[0.0f32; 3]; 3];
        for (i, row) in rgb_cam.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = libraw_get_rgb_cam(lr, i as c_int, j as c_int);
            }
        }

        let mut errc: c_int = 0;
        let img = libraw_dcraw_make_mem_image(lr, &mut errc);
        if img.is_null() {
            return Err(Error::LibRaw { code: errc, message: err_message(errc) });
        }
        let (width, height, pixels) = take_image(img)?;

        let wb = [
            libraw_get_cam_mul(lr, 0),
            libraw_get_cam_mul(lr, 1),
            libraw_get_cam_mul(lr, 2),
        ];
        let color_maximum = libraw_get_color_maximum(lr);
        // 尺寸一律经由 LibRaw 的访问接口取得，不用自己解析标签
        let iwidth = libraw_get_iwidth(lr).max(0) as usize;
        let iheight = libraw_get_iheight(lr).max(0) as usize;
        let raw_width = libraw_get_raw_width(lr).max(0) as usize;
        let raw_height = libraw_get_raw_height(lr).max(0) as usize;

        // 输出尺寸应等于 iwidth×iheight，**或者是它的转置**。
        //
        // 解码层会按 EXIF 方向把输出旋转，而 `iwidth`/`iheight` 报的是未旋转的尺寸：
        // 竖拍照片的输出是 5520×8280，而 iwidth×iheight 是 8280×5520。早前只接受
        // 相等，于是**每一张竖拍都会失败**——而竖构图在实际拍摄中很常见。
        //
        // 两种都不满足才说明解码层行为有变，那时才值得暴露。
        let same = iwidth == width && iheight == height;
        let rotated = iwidth == height && iheight == width;
        debug_assert!(
            same || rotated,
            "iwidth×iheight ({iwidth}×{iheight}) 与实际输出 ({width}×{height}) 既不相等也不互为转置"
        );

        let d = Decoded {
            width,
            height,
            rotated,
            pixels,
            wb,
            rgb_cam,
            color_maximum,
            demosaic: opts.demosaic,
            output,
            raw_width,
            raw_height,
        };
        Ok(d)
    }
}

/// 只读出有效像素尺寸，不做去马赛克与色彩转换。
///
/// 用于校验参考导出是否与原始文件同源——完整解码一张 45 MP 的 NEF 要十几秒，
/// 而尺寸在 `unpack` 之后即可取得。
pub fn read_dimensions(path: &Path) -> Result<(usize, usize)> {
    if !path.is_file() {
        return Err(Error::Io(format!("文件不存在：{}", path.display())));
    }
    let session = Session(unsafe { libraw_init(0) });
    if session.0.is_null() {
        return Err(Error::NoHandle);
    }
    let lr = session.0;
    unsafe {
        open_raw(lr, path)?;
        check(libraw_unpack(lr))?;
        Ok((
            libraw_get_raw_width(lr).max(0) as usize,
            libraw_get_raw_height(lr).max(0) as usize,
        ))
    }
}

/// 读出**完整的四元**白平衡系数 `[R, G1, B, G2]`，不做去马赛克与色彩转换。
///
/// # 为什么要看 G1 与 G2
///
/// [`read_wb`] 只返回前三项，把 `G2` 丢掉了。但尼康的白平衡**微调**（A-B 与 G-M 两轴，
/// 即 A2.0 / M2.25 这类偏移）在四元组里正是通过 **G1 与 G2 是否相等**体现的——
/// 而这是 `B/G` 一个比值**捕获不到**的一维。
///
/// 实测证据（本项目）：`DSC_0141`（K 4830）与 `DSC_0567`（K 5550 且带调色偏移）
/// 的 `B/G` 几乎相同（1.49609 vs 1.49414），而 `R/G` 差 21%（1.75000 vs 2.11133）。
/// **同一 `B/G` 对应两种不同的白平衡**，所以拿 `B/G` 当唯一自变量必然出错。
pub fn read_wb4(path: &Path) -> Result<[f32; 4]> {
    if !path.is_file() {
        return Err(Error::Io(format!("文件不存在：{}", path.display())));
    }
    let session = Session(unsafe { libraw_init(0) });
    if session.0.is_null() {
        return Err(Error::NoHandle);
    }
    let lr = session.0;
    unsafe {
        open_raw(lr, path)?;
        check(libraw_unpack(lr))?;
        Ok([
            libraw_get_cam_mul(lr, 0),
            libraw_get_cam_mul(lr, 1),
            libraw_get_cam_mul(lr, 2),
            libraw_get_cam_mul(lr, 3),
        ])
    }
}

/// 只读出白平衡系数，不做去马赛克与色彩转换。
///
/// 用于快速扫描一批文件找出白平衡不同的样本——完整解码一张 45 MP 的 NEF 要十几秒，
/// 而这一步只到 `unpack` 为止。
pub fn read_wb(path: &Path) -> Result<[f32; 3]> {
    if !path.is_file() {
        return Err(Error::Io(format!("文件不存在：{}", path.display())));
    }
    let session = Session(unsafe { libraw_init(0) });
    if session.0.is_null() {
        return Err(Error::NoHandle);
    }
    let lr = session.0;
    unsafe {
        open_raw(lr, path)?;
        check(libraw_unpack(lr))?;
        Ok([
            libraw_get_cam_mul(lr, 0),
            libraw_get_cam_mul(lr, 1),
            libraw_get_cam_mul(lr, 2),
        ])
    }
}

/// 解码一张 RAW，输出**线性相机空间**的 16 位 RGB。
///
/// 结果仍在相机空间、且不含相机内观感——转入工作空间是 [`crate::color`] 的职责。
///
/// # 为什么解码两遍
///
/// `rgb_cam` 只有在输出色彩空间**不是**相机空间时才是有效变换（相机空间下它是
/// 恒等）。而 `rgb_cam` 又在 `dcraw_process` 内部才计算。因此先跑一遍 ProPhoto
/// 取矩阵，再跑一遍取相机空间像素。这是 LibRaw 的 C API 所决定的代价。
pub fn decode_camera_linear(path: &Path, opts: &Options) -> Result<Decoded> {
    let pro = decode_with_output(path, opts, OutputColor::ProPhoto)?;
    let mut cam = decode_with_output(path, opts, OutputColor::Camera)?;
    cam.rgb_cam = pro.rgb_cam;
    Ok(cam)
}

/// 解码一张 RAW，**直接**输出工作空间（ProPhoto）线性像素。
///
/// 与「取相机空间像素 + 本管线自行施加矩阵」互为对照，用于交叉验证矩阵施加
/// 是否与解码层内部一致。日常渲染不应走这条——它会就地裁到 0..65535。
pub fn decode_working_space(path: &Path, opts: &Options) -> Result<Decoded> {
    decode_with_output(path, opts, OutputColor::ProPhoto)
}

/// 以任意输出色彩空间解码，供诊断与交叉验证使用。
pub fn decode_with_output_for_test(
    path: &Path,
    opts: &Options,
    output: OutputColor,
) -> Result<Decoded> {
    decode_with_output(path, opts, output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_against_libraw() {
        let v = version_string();
        assert_eq!(version_number() >> 16, 0, "LibRaw 主版本应为 0，得到 {v}");
        assert!(version_number() > 0, "版本号应非零");
        eprintln!("LibRaw {v} (0x{:06X})", version_number());
    }

    #[test]
    fn missing_file_is_io_error() {
        let e = decode_camera_linear(Path::new("no-such-file.nef"), &Options::default());
        assert!(matches!(e, Err(Error::Io(_))), "得到 {e:?}");
    }

    /// 损坏文件必须给出**确定的失败**，而不是崩溃或产出图像。
    ///
    /// 断言不限定具体错误变体：解码层在不同情况下可能报"不支持的格式"，也可能在
    /// 更早一步失败，两者都是合格的"确定失败"。
    ///
    /// # 这条测试曾经偶发失败，根因不在解码
    ///
    /// 失败点是 `std::fs::write` 本身，不是解码。根因有两条，都已修：
    ///
    /// 1. 临时文件名原先只带进程号——`cargo test` 并行跑多个测试二进制时可能相撞。
    ///    现掺入纳秒时间戳。
    /// 2. **本机的 `%TEMP%` 不可写**：`std::fs::write` 报 `拒绝访问 (os error 5)`。
    ///    沙箱有时把 `TEMP` 指到会话专有的可写目录、有时指到真实 temp，于是表现为
    ///    "偶发"。现改为写进 cargo 的 `target/` 目录——那是构建产物目录，必然可写，
    ///    且已在 `.gitignore` 里。
    #[test]
    fn garbage_file_fails_explicitly() {
        use std::time::{SystemTime, UNIX_EPOCH};
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);

        // 不用 temp_dir()：本机的 %TEMP% 拒绝写入（os error 5）
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
        std::fs::create_dir_all(&dir).expect("应能创建 target 目录");
        let p = dir.join(format!(
            "garbage-{}-{}.nef",
            std::process::id(),
            nanos
        ));
        if let Err(e) = std::fs::write(&p, vec![0xAAu8; 4096]) {
            panic!("写入临时文件失败：{e}（路径 {}）", p.display());
        }
        assert!(p.is_file(), "临时文件应已写入：{}", p.display());
        let e = decode_camera_linear(&p, &Options::default());
        let _ = std::fs::remove_file(&p);

        match e {
            Err(err) => {
                let msg = err.to_string();
                eprintln!("损坏文件的失败：{err:?} / {msg}");
                assert!(!msg.trim().is_empty(), "失败必须带可读的原因");
            }
            Ok(d) => panic!(
                "损坏文件不应解码成功，却得到 {}×{} 的图像",
                d.width, d.height
            ),
        }
    }

    #[test]
    fn demosaic_names_are_distinct() {
        let all = [
            Demosaic::Linear,
            Demosaic::Vng,
            Demosaic::Ppg,
            Demosaic::Ahd,
            Demosaic::Dcb,
            Demosaic::AhdMedian,
            Demosaic::Dht,
            Demosaic::Aahd,
        ];
        let mut names: Vec<&str> = all.iter().map(|d| d.name()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), all.len());
    }
}
