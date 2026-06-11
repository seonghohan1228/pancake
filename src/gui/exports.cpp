#include "exports.hpp"

#include <commdlg.h>

#include <fstream>

// Third-party header compiled with warnings-as-errors disabled locally.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-function"
#elif defined(_MSC_VER)
#pragma warning(push, 3)
#endif
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

namespace exports {

bool write_csv_columns(const std::filesystem::path& path,
                       const std::vector<std::string>& headers,
                       const std::vector<std::vector<double>>& columns,
                       std::string& error) {
    std::ofstream file(path, std::ios::trunc);
    if (!file.is_open()) {
        error = "Cannot write " + path.string();
        return false;
    }
    file.precision(12);
    for (size_t i = 0; i < headers.size(); ++i) {
        file << (i ? "," : "") << headers[i];
    }
    file << "\n";
    size_t row_count = 0;
    for (const auto& column : columns) row_count = std::max(row_count, column.size());
    for (size_t row = 0; row < row_count; ++row) {
        for (size_t i = 0; i < columns.size(); ++i) {
            if (i) file << ",";
            if (row < columns[i].size()) file << columns[i][row];
        }
        file << "\n";
    }
    return file.good();
}

bool write_csv_rows(const std::filesystem::path& path,
                    const std::vector<std::string>& headers,
                    const std::vector<std::vector<std::string>>& rows,
                    std::string& error) {
    std::ofstream file(path, std::ios::trunc);
    if (!file.is_open()) {
        error = "Cannot write " + path.string();
        return false;
    }
    const auto write_cell = [&](const std::string& cell) {
        const bool needs_quotes = cell.find_first_of(",\"\n") != std::string::npos;
        if (!needs_quotes) {
            file << cell;
            return;
        }
        file << '"';
        for (char ch : cell) {
            if (ch == '"') file << '"';
            file << ch;
        }
        file << '"';
    };
    for (size_t i = 0; i < headers.size(); ++i) {
        if (i) file << ",";
        write_cell(headers[i]);
    }
    file << "\n";
    for (const auto& row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (i) file << ",";
            write_cell(row[i]);
        }
        file << "\n";
    }
    return file.good();
}

bool write_png_rgba(const std::filesystem::path& path, int width, int height,
                    const unsigned char* pixels, std::string& error) {
    // stb takes a char* path; go through a UTF-8 conversion-safe FILE handle.
    FILE* file = _wfopen(path.c_str(), L"wb");
    if (!file) {
        error = "Cannot write " + path.string();
        return false;
    }
    const auto write_callback = [](void* context, void* data, int size) {
        fwrite(data, 1, static_cast<size_t>(size), static_cast<FILE*>(context));
    };
    const int ok =
        stbi_write_png_to_func(write_callback, file, width, height, 4, pixels, width * 4);
    fclose(file);
    if (!ok) {
        error = "PNG encode failed for " + path.string();
        return false;
    }
    return true;
}

std::filesystem::path save_file_dialog(HWND owner, const wchar_t* filter,
                                       const wchar_t* default_ext,
                                       const std::wstring& suggested_name) {
    wchar_t buffer[MAX_PATH] = {};
    suggested_name.copy(buffer, MAX_PATH - 1);
    OPENFILENAMEW dialog{};
    dialog.lStructSize = sizeof(dialog);
    dialog.hwndOwner = owner;
    dialog.lpstrFilter = filter;
    dialog.lpstrFile = buffer;
    dialog.nMaxFile = MAX_PATH;
    dialog.lpstrDefExt = default_ext;
    dialog.Flags = OFN_OVERWRITEPROMPT | OFN_NOCHANGEDIR;
    if (!GetSaveFileNameW(&dialog)) return {};
    return std::filesystem::path(buffer);
}

std::filesystem::path open_file_dialog(HWND owner, const wchar_t* filter) {
    wchar_t buffer[MAX_PATH] = {};
    OPENFILENAMEW dialog{};
    dialog.lStructSize = sizeof(dialog);
    dialog.hwndOwner = owner;
    dialog.lpstrFilter = filter;
    dialog.lpstrFile = buffer;
    dialog.nMaxFile = MAX_PATH;
    dialog.Flags = OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetOpenFileNameW(&dialog)) return {};
    return std::filesystem::path(buffer);
}

}  // namespace exports
