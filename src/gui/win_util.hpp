#pragma once

// Small Win32 helpers shared by the GUI. Header-only.

#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <filesystem>
#include <string>

namespace winutil {

inline std::wstring utf8_to_wide(const std::string& text) {
    if (text.empty()) return {};
    const int needed = MultiByteToWideChar(CP_UTF8, 0, text.data(),
                                           static_cast<int>(text.size()), nullptr, 0);
    std::wstring wide(static_cast<size_t>(needed), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, text.data(), static_cast<int>(text.size()),
                        wide.data(), needed);
    return wide;
}

inline std::string wide_to_utf8(const std::wstring& text) {
    if (text.empty()) return {};
    const int needed = WideCharToMultiByte(CP_UTF8, 0, text.data(),
                                           static_cast<int>(text.size()), nullptr, 0,
                                           nullptr, nullptr);
    std::string utf8(static_cast<size_t>(needed), '\0');
    WideCharToMultiByte(CP_UTF8, 0, text.data(), static_cast<int>(text.size()),
                        utf8.data(), needed, nullptr, nullptr);
    return utf8;
}

inline std::string path_to_utf8(const std::filesystem::path& path) {
    return wide_to_utf8(path.wstring());
}

inline std::filesystem::path utf8_to_path(const std::string& text) {
    return std::filesystem::path(utf8_to_wide(text));
}

/// Directory containing the running executable.
inline std::filesystem::path executable_dir() {
    wchar_t buffer[MAX_PATH] = {};
    GetModuleFileNameW(nullptr, buffer, MAX_PATH);
    return std::filesystem::path(buffer).parent_path();
}

}  // namespace winutil
