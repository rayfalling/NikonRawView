//! 构建 vendored 的 LibRaw。
//!
//! # 为什么是 nmake 而不是 MSBuild
//!
//! 见 `openspec/changes/add-nikon-render-pipeline/design.md` 决策二。简言之，MSBuild 有两个
//! 独立障碍：
//!
//! 1. `buildfiles/libraw.vcxproj` 钉死 Windows SDK `10.0.18362.0` 与平台工具集 `v142`
//! 2. 覆盖这两项后仍会失败于 `MSB6001`：进程环境块里同时存在大小写两份代理变量时，
//!    MSBuild 的 .NET 不区分大小写字典在构造子进程环境时直接抛异常
//!
//! nmake 是纯 Win32 工具，不走 .NET，因此不受第二条影响；第一条则完全绕开。
//! LibRaw 0.22 **没有 CMakeLists.txt**，nmake 是其官方构建路径之一。
//!
//! # 产物
//!
//! `third_party/libraw/lib/libraw_static.lib`（静态库，约 4.3 MB）。
//! 中间产物落在 `third_party/libraw/object/` 与 `bin/`，二者不入库。

use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let root = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR 未设置"));
    let libraw = root.join("third_party").join("libraw");

    println!("cargo:rerun-if-changed=build.rs");
    // LibRaw 的源码变动需要重新构建
    println!("cargo:rerun-if-changed={}", libraw.join("src").display());
    println!("cargo:rerun-if-changed={}", libraw.join("Makefile.msvc").display());

    if !libraw.join("libraw").join("libraw.h").exists() {
        panic!(
            "找不到 vendored 的 LibRaw 源码：{}\n\
             请确认 third_party/libraw/ 已随仓库检出。",
            libraw.join("libraw").display()
        );
    }

    // Windows 与非 Windows 各有一份同名实现，由 cfg 选择
    build_libraw(&libraw);

    let libdir = libraw.join("lib");
    println!("cargo:rustc-link-search=native={}", libdir.display());
    println!("cargo:rustc-link-lib=static=libraw_static");

    if cfg!(windows) {
        // LibRaw 在 Windows 上需要 Winsock
        println!("cargo:rustc-link-lib=ws2_32");
    } else {
        // C++ 运行库；MSVC 侧由 rustc 自行处理
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }
}

#[cfg(windows)]
fn build_libraw(libraw: &Path) {
    let static_lib = libraw.join("lib").join("libraw_static.lib");
    if static_lib.exists() {
        return;
    }

    // object/ 与 bin/ 不入库，首次构建需自行创建——缺了它们 nmake 会以
    // 「Cannot open compiler generated file」失败
    for d in ["object", "bin"] {
        std::fs::create_dir_all(libraw.join(d))
            .unwrap_or_else(|e| panic!("无法创建 {}: {e}", libraw.join(d).display()));
    }

    let vcvars = find_vcvars().unwrap_or_else(|| {
        panic!(
            "找不到 vcvars64.bat。\n\
             它随 Visual Studio 的「使用 C++ 的桌面开发」工作负载安装。\n\
             若已安装仍找不到，请设置环境变量 VCVARS64 指向该文件的完整路径。"
        )
    });

    // 通过批处理调用：需要先 call vcvars64.bat 建立 MSVC 环境，nmake 才能找到 cl.exe
    let bat = libraw.join("build_libraw.bat");
    let script = format!(
        "@echo off\r\n\
         set HTTP_PROXY=\r\n\
         set http_proxy=\r\n\
         set HTTPS_PROXY=\r\n\
         set https_proxy=\r\n\
         call \"{}\" >nul\r\n\
         cd /d \"{}\"\r\n\
         nmake /f Makefile.msvc lib\\libraw_static.lib\r\n\
         exit /b %ERRORLEVEL%\r\n",
        vcvars.display(),
        libraw.display()
    );
    std::fs::write(&bat, script).expect("无法写入构建批处理");

    let out = Command::new("cmd")
        .arg("/c")
        .arg(&bat)
        .output()
        .expect("无法启动 cmd.exe");

    let _ = std::fs::remove_file(&bat);

    if !out.status.success() || !static_lib.exists() {
        let stdout = String::from_utf8_lossy(&out.stdout);
        let tail: Vec<&str> = stdout.lines().filter(|l| !l.trim().is_empty()).rev().take(25).collect();
        panic!(
            "LibRaw 构建失败（退出码 {:?}）。\n\
             构建目录：{}\n\
             输出末尾（倒序）：\n{}",
            out.status.code(),
            libraw.display(),
            tail.iter().rev().map(|l| format!("  {l}")).collect::<Vec<_>>().join("\n")
        );
    }
}

#[cfg(windows)]
fn find_vcvars() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("VCVARS64") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }

    // 首选 vswhere：它能处理多版本共存与自定义安装路径
    let vswhere = PathBuf::from(
        std::env::var("ProgramFiles(x86)").unwrap_or_else(|_| r"C:\Program Files (x86)".into()),
    )
    .join("Microsoft Visual Studio")
    .join("Installer")
    .join("vswhere.exe");

    if vswhere.is_file() {
        let out = Command::new(&vswhere)
            .args([
                "-latest",
                "-prerelease",
                "-requires",
                "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                "-find",
                r"VC\Auxiliary\Build\vcvars64.bat",
            ])
            .output()
            .ok()?;
        if out.status.success() {
            let s = String::from_utf8_lossy(&out.stdout);
            if let Some(line) = s.lines().map(str::trim).find(|l| !l.is_empty()) {
                let p = PathBuf::from(line);
                if p.is_file() {
                    return Some(p);
                }
            }
        }
    }

    // 退路：扫描常见安装位置
    for base in [
        r"C:\Program Files\Microsoft Visual Studio",
        r"C:\Program Files (x86)\Microsoft Visual Studio",
    ] {
        let base = PathBuf::from(base);
        if !base.is_dir() {
            continue;
        }
        for year in ["2022", "2019", "18", "17"] {
            for ed in ["Enterprise", "Professional", "Community", "BuildTools"] {
                let p = base
                    .join(year)
                    .join(ed)
                    .join(r"VC\Auxiliary\Build\vcvars64.bat");
                if p.is_file() {
                    return Some(p);
                }
            }
        }
    }
    None
}

#[cfg(not(windows))]
fn build_libraw(libraw: &Path) {
    let static_lib = libraw.join("lib").join("libraw.a");
    if static_lib.exists() {
        return;
    }
    std::fs::create_dir_all(libraw.join("object")).ok();
    let jobs = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    let out = Command::new("make")
        .current_dir(libraw)
        .args(["-f", "Makefile.dist", &format!("-j{jobs}"), "lib/libraw.a"])
        .output()
        .expect("无法启动 make");
    if !out.status.success() || !static_lib.exists() {
        panic!(
            "LibRaw 构建失败（非 Windows）：\n{}",
            String::from_utf8_lossy(&out.stdout)
        );
    }
}
