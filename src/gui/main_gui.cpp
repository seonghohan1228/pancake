// Platform layer: Win32 window, Direct3D 11 device/swapchain, Dear ImGui
// (docking) + ImPlot lifecycle, per-monitor-v2 DPI handling, frame loop and
// plot-to-PNG capture. Application logic lives in GuiApp.

#include "gui_app.hpp"

#include "exports.hpp"
#include "imgui.h"
#include "imgui_impl_dx11.h"
#include "imgui_impl_win32.h"
#include "implot.h"
#include "settings.hpp"
#include "theme.hpp"
#include "win_util.hpp"

#include <d3d11.h>

#include <algorithm>
#include <string>
#include <vector>

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hwnd, UINT msg, WPARAM wparam,
                                                             LPARAM lparam);

namespace {

struct Platform {
    HWND hwnd = nullptr;
    ID3D11Device* device = nullptr;
    ID3D11DeviceContext* context = nullptr;
    IDXGISwapChain* swap_chain = nullptr;
    ID3D11RenderTargetView* render_target = nullptr;
    UINT resize_width = 0;
    UINT resize_height = 0;
    float dpi_scale = 1.0f;
    bool dpi_changed = false;
    bool occluded = false;
    GuiApp* app = nullptr;
};

void create_render_target(Platform& platform) {
    ID3D11Texture2D* back_buffer = nullptr;
    platform.swap_chain->GetBuffer(0, IID_PPV_ARGS(&back_buffer));
    if (back_buffer) {
        platform.device->CreateRenderTargetView(back_buffer, nullptr, &platform.render_target);
        back_buffer->Release();
    }
}

void release_render_target(Platform& platform) {
    if (platform.render_target) {
        platform.render_target->Release();
        platform.render_target = nullptr;
    }
}

bool create_device(Platform& platform) {
    DXGI_SWAP_CHAIN_DESC swap_desc{};
    swap_desc.BufferCount = 2;
    swap_desc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    swap_desc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    swap_desc.OutputWindow = platform.hwnd;
    swap_desc.SampleDesc.Count = 1;
    swap_desc.Windowed = TRUE;
    swap_desc.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
    swap_desc.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;

    const D3D_FEATURE_LEVEL levels[] = {D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0};
    D3D_FEATURE_LEVEL created_level;
    HRESULT result = D3D11CreateDeviceAndSwapChain(
        nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 0, levels, 2, D3D11_SDK_VERSION, &swap_desc,
        &platform.swap_chain, &platform.device, &created_level, &platform.context);
    if (result == DXGI_ERROR_UNSUPPORTED) {
        result = D3D11CreateDeviceAndSwapChain(
            nullptr, D3D_DRIVER_TYPE_WARP, nullptr, 0, levels, 2, D3D11_SDK_VERSION, &swap_desc,
            &platform.swap_chain, &platform.device, &created_level, &platform.context);
    }
    if (FAILED(result)) return false;
    create_render_target(platform);
    return true;
}

void destroy_device(Platform& platform) {
    release_render_target(platform);
    if (platform.swap_chain) platform.swap_chain->Release();
    if (platform.context) platform.context->Release();
    if (platform.device) platform.device->Release();
    platform.swap_chain = nullptr;
    platform.context = nullptr;
    platform.device = nullptr;
}

/// Copies the backbuffer region [min, max) to a PNG file.
void capture_region_to_png(Platform& platform, const GuiApp::CaptureRequest& request) {
    ID3D11Texture2D* back_buffer = nullptr;
    platform.swap_chain->GetBuffer(0, IID_PPV_ARGS(&back_buffer));
    if (!back_buffer) {
        platform.app->capture_finished(false, "no backbuffer");
        return;
    }
    D3D11_TEXTURE2D_DESC desc{};
    back_buffer->GetDesc(&desc);
    desc.Usage = D3D11_USAGE_STAGING;
    desc.BindFlags = 0;
    desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
    desc.MiscFlags = 0;
    ID3D11Texture2D* staging = nullptr;
    platform.device->CreateTexture2D(&desc, nullptr, &staging);
    if (!staging) {
        back_buffer->Release();
        platform.app->capture_finished(false, "staging texture failed");
        return;
    }
    platform.context->CopyResource(staging, back_buffer);
    back_buffer->Release();

    D3D11_MAPPED_SUBRESOURCE mapped{};
    if (FAILED(platform.context->Map(staging, 0, D3D11_MAP_READ, 0, &mapped))) {
        staging->Release();
        platform.app->capture_finished(false, "map failed");
        return;
    }

    const int x0 = std::clamp(static_cast<int>(request.min_x) - 4, 0, static_cast<int>(desc.Width));
    const int y0 = std::clamp(static_cast<int>(request.min_y) - 4, 0, static_cast<int>(desc.Height));
    const int x1 = std::clamp(static_cast<int>(request.max_x) + 4, 0, static_cast<int>(desc.Width));
    const int y1 = std::clamp(static_cast<int>(request.max_y) + 4, 0, static_cast<int>(desc.Height));
    const int width = x1 - x0;
    const int height = y1 - y0;
    if (width <= 0 || height <= 0) {
        platform.context->Unmap(staging, 0);
        staging->Release();
        platform.app->capture_finished(false, "empty capture region");
        return;
    }

    std::vector<unsigned char> pixels(static_cast<size_t>(width) * height * 4);
    const auto* source = static_cast<const unsigned char*>(mapped.pData);
    for (int row = 0; row < height; ++row) {
        const unsigned char* line = source + static_cast<size_t>(y0 + row) * mapped.RowPitch +
                                    static_cast<size_t>(x0) * 4;
        unsigned char* target = pixels.data() + static_cast<size_t>(row) * width * 4;
        for (int col = 0; col < width; ++col) {
            target[col * 4 + 0] = line[col * 4 + 0];  // backbuffer is RGBA8
            target[col * 4 + 1] = line[col * 4 + 1];
            target[col * 4 + 2] = line[col * 4 + 2];
            target[col * 4 + 3] = 255;  // exported image is opaque
        }
    }
    platform.context->Unmap(staging, 0);
    staging->Release();

    std::string error;
    if (exports::write_png_rgba(request.path, width, height, pixels.data(), error)) {
        platform.app->capture_finished(true, winutil::path_to_utf8(request.path));
    } else {
        platform.app->capture_finished(false, error);
    }
}

LRESULT WINAPI wnd_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam) {
    if (ImGui_ImplWin32_WndProcHandler(hwnd, msg, wparam, lparam)) return 1;
    auto* platform = reinterpret_cast<Platform*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    switch (msg) {
        case WM_SIZE:
            if (platform && wparam != SIZE_MINIMIZED) {
                platform->resize_width = LOWORD(lparam);
                platform->resize_height = HIWORD(lparam);
            }
            return 0;
        case WM_DPICHANGED: {
            if (platform) {
                platform->dpi_scale = static_cast<float>(HIWORD(wparam)) / 96.0f;
                platform->dpi_changed = true;
            }
            const RECT* suggested = reinterpret_cast<RECT*>(lparam);
            SetWindowPos(hwnd, nullptr, suggested->left, suggested->top,
                         suggested->right - suggested->left, suggested->bottom - suggested->top,
                         SWP_NOZORDER | SWP_NOACTIVATE);
            return 0;
        }
        case WM_GETMINMAXINFO: {
            auto* info = reinterpret_cast<MINMAXINFO*>(lparam);
            const float scale = platform ? platform->dpi_scale : 1.0f;
            info->ptMinTrackSize.x = static_cast<LONG>(880 * scale);
            info->ptMinTrackSize.y = static_cast<LONG>(560 * scale);
            return 0;
        }
        case WM_CLOSE:
            if (platform && platform->app && platform->app->solver_running()) {
                if (MessageBoxW(hwnd,
                                L"A solve is still running. Stop it and exit?",
                                L"Pancake GUI", MB_YESNO | MB_ICONWARNING) != IDYES) {
                    return 0;
                }
            }
            DestroyWindow(hwnd);
            return 0;
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
        default: break;
    }
    return DefWindowProcW(hwnd, msg, wparam, lparam);
}

}  // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show_command) {
    ImGui_ImplWin32_EnableDpiAwareness();  // manifest declares PMv2; this is the fallback

    Platform platform;
    GuiApp app;
    platform.app = &app;
    app.paths = AppPaths::resolve();
    app.settings.load(app.paths.settings_file);

    WNDCLASSEXW window_class{};
    window_class.cbSize = sizeof(window_class);
    window_class.style = CS_HREDRAW | CS_VREDRAW;
    window_class.lpfnWndProc = wnd_proc;
    window_class.hInstance = instance;
    window_class.hIcon = LoadIconW(instance, L"IDI_PANCAKE");
    window_class.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    window_class.lpszClassName = L"PancakeGuiMain";
    RegisterClassExW(&window_class);

    platform.hwnd = CreateWindowExW(
        0, window_class.lpszClassName, L"Pancake - journal bearing solver",
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, app.settings.window_width,
        app.settings.window_height, nullptr, nullptr, instance, nullptr);
    SetWindowLongPtrW(platform.hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(&platform));
    platform.dpi_scale = ImGui_ImplWin32_GetDpiScaleForHwnd(platform.hwnd);

    if (!create_device(platform)) {
        MessageBoxW(platform.hwnd,
                    L"Direct3D 11 initialization failed (no hardware or WARP device).",
                    L"Pancake GUI", MB_OK | MB_ICONERROR);
        DestroyWindow(platform.hwnd);
        UnregisterClassW(window_class.lpszClassName, instance);
        return 1;
    }

    ShowWindow(platform.hwnd,
               app.settings.window_maximized ? SW_SHOWMAXIMIZED : show_command);
    UpdateWindow(platform.hwnd);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard | ImGuiConfigFlags_DockingEnable;
    static std::string ini_path_utf8 = winutil::path_to_utf8(app.paths.imgui_ini);
    io.IniFilename = ini_path_utf8.c_str();

    ImGui_ImplWin32_Init(platform.hwnd);
    ImGui_ImplDX11_Init(platform.device, platform.context);

    theme::apply_style(platform.dpi_scale);
    app.fonts = theme::build_fonts(platform.dpi_scale);
    app.hwnd = platform.hwnd;
    app.init();

    bool done = false;
    while (!done) {
        MSG msg;
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
            if (msg.message == WM_QUIT) done = true;
        }
        if (done || app.wants_exit) break;

        if (platform.occluded &&
            platform.swap_chain->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) {
            Sleep(10);
            continue;
        }
        platform.occluded = false;

        if (platform.resize_width != 0 && platform.resize_height != 0) {
            release_render_target(platform);
            platform.swap_chain->ResizeBuffers(0, platform.resize_width, platform.resize_height,
                                               DXGI_FORMAT_UNKNOWN, 0);
            platform.resize_width = platform.resize_height = 0;
            create_render_target(platform);
        }

        if (platform.dpi_changed) {
            platform.dpi_changed = false;
            ImGui_ImplDX11_InvalidateDeviceObjects();
            theme::apply_style(platform.dpi_scale);
            app.fonts = theme::build_fonts(platform.dpi_scale);
            ImGui_ImplDX11_CreateDeviceObjects();
        }

        ImGui_ImplDX11_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();
        app.draw();
        ImGui::Render();

        constexpr float kClear[4] = {0.08f, 0.085f, 0.095f, 1.0f};
        platform.context->OMSetRenderTargets(1, &platform.render_target, nullptr);
        platform.context->ClearRenderTargetView(platform.render_target, kClear);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

        if (app.pending_capture) {
            const GuiApp::CaptureRequest request = *app.pending_capture;
            app.pending_capture.reset();
            capture_region_to_png(platform, request);
        }

        const HRESULT present_result = platform.swap_chain->Present(1, 0);
        platform.occluded = present_result == DXGI_STATUS_OCCLUDED;
    }

    // Persist window placement for the next session.
    WINDOWPLACEMENT placement{};
    placement.length = sizeof(placement);
    if (GetWindowPlacement(platform.hwnd, &placement)) {
        app.settings.window_maximized = placement.showCmd == SW_SHOWMAXIMIZED;
        app.settings.window_width =
            std::max(640L, placement.rcNormalPosition.right - placement.rcNormalPosition.left);
        app.settings.window_height =
            std::max(480L, placement.rcNormalPosition.bottom - placement.rcNormalPosition.top);
    }
    app.shutdown();

    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    destroy_device(platform);
    if (IsWindow(platform.hwnd)) DestroyWindow(platform.hwnd);
    UnregisterClassW(window_class.lpszClassName, instance);
    return 0;
}
