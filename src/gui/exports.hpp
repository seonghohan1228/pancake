#pragma once

// CSV/PNG export helpers and Win32 file dialogs.

#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <filesystem>
#include <string>
#include <vector>

namespace exports {

/// Numeric column-oriented CSV (curves, fields). Columns may differ in length;
/// short columns leave trailing cells empty.
bool write_csv_columns(const std::filesystem::path& path,
                       const std::vector<std::string>& headers,
                       const std::vector<std::vector<double>>& columns,
                       std::string& error);

/// Row-oriented CSV of preformatted cells (tables).
bool write_csv_rows(const std::filesystem::path& path,
                    const std::vector<std::string>& headers,
                    const std::vector<std::vector<std::string>>& rows,
                    std::string& error);

/// RGBA8 image to PNG via stb_image_write.
bool write_png_rgba(const std::filesystem::path& path, int width, int height,
                    const unsigned char* pixels, std::string& error);

/// Returns an empty path if the user cancels.
std::filesystem::path save_file_dialog(HWND owner, const wchar_t* filter,
                                       const wchar_t* default_ext,
                                       const std::wstring& suggested_name);
std::filesystem::path open_file_dialog(HWND owner, const wchar_t* filter);

}  // namespace exports
