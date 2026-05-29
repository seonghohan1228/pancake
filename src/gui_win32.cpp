#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <windows.h>
#include <commctrl.h>
#include <richedit.h>
#include <shellapi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cwctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.hpp"
#include "resource.h"

#ifdef _MSC_VER
#pragma comment(linker, "\"/manifestdependency:type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")
#endif

namespace {
namespace fs = std::filesystem;

constexpr wchar_t kWindowClass[] = L"PancakeGuiWindow";
constexpr wchar_t kPreviewClass[] = L"PancakePreviewCanvas";
constexpr wchar_t kSplitterClass[] = L"PancakeSplitter";
constexpr wchar_t kParamPanelClass[] = L"PancakeParamPanel";
constexpr UINT_PTR kProcessTimer = 1;
constexpr int kSplitterWidth = 8;
constexpr int kPaneGap = 10;
constexpr double kPi = 3.14159265358979323846264338327950288;

// Modern light palette shared by the window background and the WM_CTLCOLOR*
// handlers so labels, group boxes, and check boxes blend into one clean
// surface instead of the default system grey.
constexpr COLORREF kColPage = RGB(244, 247, 251);   // window / static / groupbox surface
constexpr COLORREF kColField = RGB(255, 255, 255);  // edit and list backgrounds
constexpr COLORREF kColText = RGB(28, 36, 46);      // primary text
constexpr COLORREF kColMuted = RGB(108, 119, 133);  // secondary text
constexpr COLORREF kColBorder = RGB(208, 215, 224); // hairline separators
constexpr COLORREF kColAccent = RGB(38, 110, 196);  // headings / highlights

enum TabIndex {
    kWorkspaceTab = 0,
    kOutputTab = 1,
    kRawTab = 2,
    kTabCount = 3,
    kParametersTab = kWorkspaceTab,
    kRunTab = kWorkspaceTab
};

enum ControlId {
    kTabsId = 1000,
    kSaveId,
    kReloadId,
    kDefaultsId,
    kRunId,
    kStopId,
    kOpenResultsId,
    kStatusId,

    kEditR = 2000,
    kEditC,
    kEditE,
    kEditL,
    kEditAttitude,
    kEditLoadAngle,
    kEditTiltX,
    kEditTiltY,
    kEditNTheta,
    kEditNZ,
    kEditEndT,
    kEditDt,
    kEditWriteInterval,
    kEditOmega,
    kEditMu,
    kEditRho,
    kEditPCav,
    kEditBulkModulus,
    kEditMaxOuter,
    kEditOuterTol,
    kEditThetaMin,
    kEditBcSouthVal,
    kEditBcSouthTheta,
    kEditBcNorthVal,
    kEditBcNorthTheta,
    kEditInlets,
    kEditInletTheta,
    kEditInletZ,
    kEditInletSize,
    kEditInletPressure,
    kEditOutputDir,
    kEditFilenamePrefix,
    kEditSolverPath,
    kEditMpiRanks,
    kEditLog,
    kEditRawConfig,

    kComboCavitation = 3000,
    kComboBcSouth,
    kComboBcNorth,
    kComboPreviewField,
    kComboPreviewStep,
    kComboInletType,
    kComboUnitR,
    kComboUnitC,
    kComboUnitE,
    kComboUnitL,
    kComboUnitAttitude,
    kComboUnitLoadAngle,
    kComboUnitTiltX,
    kComboUnitTiltY,
    kComboUnitEndT,
    kComboUnitDt,
    kComboUnitWriteInterval,
    kComboUnitOmega,
    kComboUnitInletTheta,
    kComboUnitInletZ,
    kComboUnitInletSize,
    kComboUnitInletPressure,
    kComboUnitNTheta,
    kComboUnitNZ,
    kComboUnitMu,
    kComboUnitRho,
    kComboUnitPCav,
    kComboUnitBulkModulus,
    kComboUnitBcSouthVal,
    kComboUnitBcSouthTheta,
    kComboUnitBcNorthVal,
    kComboUnitBcNorthTheta,
    kComboUnitMaxOuter,
    kComboUnitOuterTol,
    kComboUnitThetaMin,

    kCheckLogOuter = 4000,
    kCheckOutput3d,
    kCheckOutputFlat,
    kCheckFieldPressure,
    kCheckFieldFilmContent,
    kCheckFieldH,
    kCheckFieldRho,
    kCheckFieldInlet,
    kCheckFieldVelocity,
    kCheckFieldLoadX,
    kCheckFieldLoadY,
    kCheckFieldLoadZ,
    kCheckFieldTorque,

    kRefreshPreviewId = 5000,
    kOpenParaviewId,
    kPrevStepId,
    kNextStepId,
    kSyncRawId,
    kApplyRawId,
    kPreviewCanvasId,
    kRunPanelId,
    kStopPanelId
};

struct ControlRef {
    int id = 0;
    HWND hwnd = nullptr;
};

struct OutputFieldControl {
    int id = 0;
    const char* key = nullptr;
    const wchar_t* label = nullptr;
};

enum class ParamVisibility {
    Always,
    CircularInletOnly
};

enum class ParamItemKind {
    Group,
    EditRow,
    ComboRow,
    CheckRow
};

struct ParamItem {
    ParamItemKind kind = ParamItemKind::EditRow;
    ParamVisibility visibility = ParamVisibility::Always;
    HWND group = nullptr;
    HWND label = nullptr;
    HWND value = nullptr;
    HWND unit = nullptr;
};

struct PreviewData {
    bool ok = false;
    int nx = 0;
    int ny = 0;
    int step = -1;
    double min_value = 0.0;
    double max_value = 0.0;
    std::string field;
    std::wstring message = L"No result loaded.";
    std::wstring source_path;
    std::vector<double> values;
};

struct PreviewStep {
    int step = -1;
    fs::path path;
    fs::file_time_type write_time{};
};

struct SolverProcess {
    PROCESS_INFORMATION process_info{};
    HANDLE stdout_read = nullptr;
    bool running = false;
};

const std::array<OutputFieldControl, 10> kOutputFields = {{
    {kCheckFieldPressure, "pressure", L"pressure"},
    {kCheckFieldFilmContent, "film_content", L"film_content"},
    {kCheckFieldH, "h", L"h"},
    {kCheckFieldRho, "rho", L"rho"},
    {kCheckFieldInlet, "inlet_indicator", L"inlet_indicator"},
    {kCheckFieldVelocity, "velocity", L"velocity"},
    {kCheckFieldLoadX, "load_x", L"load_x"},
    {kCheckFieldLoadY, "load_y", L"load_y"},
    {kCheckFieldLoadZ, "load_z", L"load_z"},
    {kCheckFieldTorque, "friction_torque", L"friction_torque"},
}};

HWND g_hwnd = nullptr;
HWND g_tabs = nullptr;
HWND g_status = nullptr;
HWND g_run_button = nullptr;
HWND g_stop_button = nullptr;
HWND g_run_panel_button = nullptr;
HWND g_stop_panel_button = nullptr;
HWND g_preview_canvas = nullptr;
HWND g_run_result_group = nullptr;
HWND g_run_solver_group = nullptr;
HWND g_run_log_group = nullptr;
HWND g_run_splitter = nullptr;
HWND g_run_field_label = nullptr;
HWND g_run_step_label = nullptr;
HWND g_run_solver_label = nullptr;
HWND g_run_ranks_label = nullptr;
HWND g_inlet_z_label = nullptr;
HWND g_inlet_size_label = nullptr;
HWND g_param_panel = nullptr;
HWND g_parent_override = nullptr;
HFONT g_font = nullptr;
HFONT g_heading_font = nullptr;
HFONT g_mono_font = nullptr;
HBRUSH g_brush_page = nullptr;
HBRUSH g_brush_field = nullptr;
HMODULE g_richedit_module = nullptr;
bool g_richedit_available = false;
fs::path g_executable_dir;
fs::path g_config_path;
SimulationConfig g_config;
std::map<std::string, std::string> g_config_value_text;
std::string g_config_inlets_text;
SolverProcess g_solver;
PreviewData g_preview;
std::vector<PreviewStep> g_preview_steps;
bool g_dirty = false;
bool g_loading = false;
int g_active_tab = kParametersTab;
int g_run_bottom_height = 190;
int g_param_scroll_y = 0;
int g_param_content_height = 0;
int g_param_panel_width = 360;
bool g_run_splitter_dragging = false;
std::vector<ControlRef> g_control_refs;
std::array<std::vector<HWND>, kTabCount> g_tab_controls;
std::vector<ParamItem> g_param_items;

std::wstring utf8_to_wide(const std::string& text) {
    if (text.empty()) return {};

    int length = MultiByteToWideChar(
        CP_UTF8, MB_ERR_INVALID_CHARS, text.data(), static_cast<int>(text.size()), nullptr, 0);
    UINT code_page = CP_UTF8;
    DWORD flags = MB_ERR_INVALID_CHARS;

    if (length == 0) {
        code_page = CP_ACP;
        flags = 0;
        length = MultiByteToWideChar(
            code_page, flags, text.data(), static_cast<int>(text.size()), nullptr, 0);
    }

    std::wstring wide(static_cast<size_t>(std::max(length, 0)), L'\0');
    if (length > 0) {
        MultiByteToWideChar(
            code_page, flags, text.data(), static_cast<int>(text.size()), wide.data(), length);
    }
    return wide;
}

std::string wide_to_utf8(const std::wstring& text) {
    if (text.empty()) return {};

    const int length = WideCharToMultiByte(
        CP_UTF8, 0, text.data(), static_cast<int>(text.size()), nullptr, 0, nullptr, nullptr);
    std::string utf8(static_cast<size_t>(std::max(length, 0)), '\0');
    if (length > 0) {
        WideCharToMultiByte(
            CP_UTF8, 0, text.data(), static_cast<int>(text.size()), utf8.data(), length, nullptr, nullptr);
    }
    return utf8;
}

fs::path executable_directory() {
    std::vector<wchar_t> buffer(MAX_PATH);
    while (true) {
        const DWORD length = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (length == 0) return fs::current_path();
        if (length < buffer.size()) return fs::path(buffer.data()).parent_path();
        buffer.resize(buffer.size() * 2);
    }
}

std::string read_file_text(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return {};
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}

std::wstring get_window_text(HWND hwnd) {
    const int length = GetWindowTextLengthW(hwnd);
    std::vector<wchar_t> buffer(static_cast<size_t>(length) + 1);
    GetWindowTextW(hwnd, buffer.data(), static_cast<int>(buffer.size()));
    return std::wstring(buffer.data(), static_cast<size_t>(length));
}

std::wstring normalize_newlines_for_edit(const std::wstring& text) {
    std::wstring normalized;
    normalized.reserve(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] == L'\r') {
            normalized += L'\r';
            if (i + 1 < text.size() && text[i + 1] == L'\n') {
                normalized += L'\n';
                ++i;
            } else {
                normalized += L'\n';
            }
        } else if (text[i] == L'\n') {
            normalized += L"\r\n";
        } else {
            normalized += text[i];
        }
    }
    return normalized;
}

std::string get_edit_utf8(int id);

std::string trim_copy(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string normalise_key(std::string text) {
    text = trim_copy(text);
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return text;
}

std::string compact_number(double value) {
    if (value == 0.0) return "0";
    std::ostringstream out;
    out << std::setprecision(12) << std::defaultfloat << value;
    return out.str();
}

std::map<std::string, std::string> parse_config_value_text(const std::string& text) {
    std::map<std::string, std::string> values;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        const size_t comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);

        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = trim_copy(line.substr(0, eq));
        std::string value = trim_copy(line.substr(eq + 1));
        if (!key.empty()) values[key] = value;
    }
    return values;
}

std::string parse_config_inlets_text(const std::string& text) {
    std::ostringstream out;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        const size_t comment = line.find('#');
        std::string active = (comment == std::string::npos) ? line : line.substr(0, comment);
        const size_t eq = active.find('=');
        if (eq == std::string::npos) continue;

        const std::string key = normalise_key(active.substr(0, eq));
        if (key == "inlet_circular" || key == "inlet_groove") {
            out << trim_copy(active) << "\r\n";
        }
    }
    return out.str();
}

HWND find_control(int id) {
    for (const auto& ref : g_control_refs) {
        if (ref.id == id) return ref.hwnd;
    }
    return nullptr;
}

void remember_control(int id, HWND hwnd) {
    if (id != 0) g_control_refs.push_back({id, hwnd});
}

void add_to_tab(int tab, HWND hwnd) {
    if (tab >= 0 && tab < kTabCount) g_tab_controls[tab].push_back(hwnd);
}

void set_status(const std::wstring& text) {
    if (g_status != nullptr) SetWindowTextW(g_status, text.c_str());
}

void move_control(HWND hwnd, int x, int y, int width, int height, BOOL repaint = TRUE) {
    if (hwnd != nullptr) {
        MoveWindow(hwnd, x, y, std::max(width, 1), std::max(height, 1), repaint);
    }
}

void move_control_id(int id, int x, int y, int width, int height, BOOL repaint = TRUE) {
    move_control(find_control(id), x, y, width, height, repaint);
}

void bring_control_to_front(HWND hwnd) {
    if (hwnd != nullptr) {
        SetWindowPos(hwnd, HWND_TOP, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
    }
}

void bring_control_id_to_front(int id) {
    bring_control_to_front(find_control(id));
}

void raise_global_header_controls() {
    bring_control_to_front(g_run_solver_group);
    bring_control_to_front(g_run_solver_label);
    bring_control_id_to_front(kEditSolverPath);
    bring_control_to_front(g_run_ranks_label);
    bring_control_id_to_front(kEditMpiRanks);
    bring_control_to_front(g_run_button);
    bring_control_to_front(g_stop_button);
}

void raise_tab_page_controls(int tab) {
    if (g_tabs != nullptr) {
        SetWindowPos(g_tabs, HWND_BOTTOM, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
    }
    if (tab < 0 || tab >= kTabCount) return;

    // Controls are stored in creation order. Group boxes are created first
    // inside each section, so walking forward leaves labels, edits, combos,
    // buttons, and canvases above their surrounding frames and above the tab.
    for (HWND child : g_tab_controls[tab]) {
        if (child != nullptr && IsWindowVisible(child)) {
            bring_control_to_front(child);
        }
    }
    raise_global_header_controls();
}

void mark_dirty() {
    if (!g_loading) {
        g_dirty = true;
        set_status(L"Unsaved changes.");
    }
}

HWND create_window_control(
    int tab,
    const wchar_t* class_name,
    const wchar_t* text,
    DWORD style,
    DWORD ex_style,
    int x,
    int y,
    int width,
    int height,
    int id) {
    HWND parent = (g_parent_override != nullptr) ? g_parent_override : g_hwnd;
    HWND hwnd = CreateWindowExW(
        ex_style,
        class_name,
        text,
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | style,
        x,
        y,
        width,
        height,
        parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr),
        nullptr);
    SendMessageW(hwnd, WM_SETFONT, reinterpret_cast<WPARAM>(g_font), TRUE);
    remember_control(id, hwnd);
    add_to_tab(tab, hwnd);
    return hwnd;
}

HWND create_label(int tab, const wchar_t* text, int x, int y, int width, int height = 22) {
    return create_window_control(tab, L"STATIC", text, SS_LEFT, 0, x, y, width, height, 0);
}

HWND create_group(int tab, const wchar_t* text, int x, int y, int width, int height) {
    HWND group = create_window_control(tab, L"BUTTON", text, BS_GROUPBOX, 0, x, y, width, height, 0);
    SendMessageW(group, WM_SETFONT, reinterpret_cast<WPARAM>(g_heading_font), TRUE);
    return group;
}

HWND create_edit(int tab, int id, int x, int y, int width, int height = 24, DWORD extra_style = ES_AUTOHSCROLL) {
    return create_window_control(
        tab,
        L"EDIT",
        L"",
        WS_TABSTOP | extra_style,
        WS_EX_CLIENTEDGE,
        x,
        y,
        width,
        height,
        id);
}

HWND create_button(int id, const wchar_t* text, int x, int y, int width, int height = 30) {
    HWND button = CreateWindowExW(
        0,
        L"BUTTON",
        text,
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_TABSTOP | BS_PUSHBUTTON,
        x,
        y,
        width,
        height,
        g_hwnd,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr),
        nullptr);
    SendMessageW(button, WM_SETFONT, reinterpret_cast<WPARAM>(g_font), TRUE);
    remember_control(id, button);
    return button;
}

HWND create_tab_button(int tab, int id, const wchar_t* text, int x, int y, int width, int height = 30) {
    return create_window_control(
        tab,
        L"BUTTON",
        text,
        WS_TABSTOP | BS_PUSHBUTTON,
        0,
        x,
        y,
        width,
        height,
        id);
}

HWND create_checkbox(int tab, int id, const wchar_t* text, int x, int y, int width, int height = 24) {
    return create_window_control(
        tab,
        L"BUTTON",
        text,
        WS_TABSTOP | BS_AUTOCHECKBOX,
        0,
        x,
        y,
        width,
        height,
        id);
}

HWND create_combo(int tab, int id, int x, int y, int width, const std::vector<const wchar_t*>& items) {
    HWND combo = create_window_control(
        tab,
        L"COMBOBOX",
        L"",
        WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST | CBS_HASSTRINGS,
        0,
        x,
        y,
        width,
        220,
        id);
    for (const wchar_t* item : items) {
        SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(item));
    }
    if (!items.empty()) SendMessageW(combo, CB_SETCURSEL, 0, 0);
    bring_control_to_front(combo);
    InvalidateRect(combo, nullptr, TRUE);
    return combo;
}

void create_row(int tab, const wchar_t* label, int id, int x, int y, int label_width = 138, int edit_width = 128) {
    create_label(tab, label, x, y + 4, label_width);
    create_edit(tab, id, x + label_width, y, edit_width);
}

void create_unit_row(
    int tab,
    const wchar_t* label,
    int edit_id,
    int unit_id,
    int x,
    int y,
    const std::vector<const wchar_t*>& units,
    int label_width = 118,
    int edit_width = 84,
    int unit_width = 72) {
    create_label(tab, label, x, y + 4, label_width);
    create_edit(tab, edit_id, x + label_width, y, edit_width);
    create_combo(tab, unit_id, x + label_width + edit_width + 6, y - 1, unit_width, units);
}

void param_group(const wchar_t* title) {
    ParamItem item{};
    item.kind = ParamItemKind::Group;
    item.group = create_group(kWorkspaceTab, title, 0, 0, 1, 1);
    g_param_items.push_back(item);
}

void param_unit_row(
    const wchar_t* label,
    int edit_id,
    int unit_id,
    const std::vector<const wchar_t*>& units,
    ParamVisibility visibility = ParamVisibility::Always) {
    ParamItem item{};
    item.kind = ParamItemKind::EditRow;
    item.visibility = visibility;
    item.label = create_label(kWorkspaceTab, label, 0, 0, 1);
    item.value = create_edit(kWorkspaceTab, edit_id, 0, 0, 1);
    item.unit = create_combo(kWorkspaceTab, unit_id, 0, 0, 1, units);
    g_param_items.push_back(item);
}

void param_combo_row(
    const wchar_t* label,
    int combo_id,
    const std::vector<const wchar_t*>& items,
    ParamVisibility visibility = ParamVisibility::Always) {
    ParamItem item{};
    item.kind = ParamItemKind::ComboRow;
    item.visibility = visibility;
    item.label = create_label(kWorkspaceTab, label, 0, 0, 1);
    item.value = create_combo(kWorkspaceTab, combo_id, 0, 0, 1, items);
    g_param_items.push_back(item);
}

void param_check_row(int check_id, const wchar_t* text) {
    ParamItem item{};
    item.kind = ParamItemKind::CheckRow;
    item.value = create_checkbox(kWorkspaceTab, check_id, text, 0, 0, 1);
    g_param_items.push_back(item);
}

void set_edit_text(int id, const std::string& text) {
    HWND hwnd = find_control(id);
    if (hwnd != nullptr) SetWindowTextW(hwnd, normalize_newlines_for_edit(utf8_to_wide(text)).c_str());
}

void set_edit_text(int id, double value) {
    set_edit_text(id, compact_number(value));
}

void set_edit_text(int id, int value) {
    set_edit_text(id, std::to_string(value));
}

void set_edit_from_config_text(int id, const char* key, double fallback) {
    auto it = g_config_value_text.find(key);
    set_edit_text(id, it != g_config_value_text.end() ? it->second : compact_number(fallback));
}

void set_edit_from_config_text(int id, const char* key, int fallback) {
    auto it = g_config_value_text.find(key);
    set_edit_text(id, it != g_config_value_text.end() ? it->second : std::to_string(fallback));
}

bool combo_has_selection(int id) {
    HWND hwnd = find_control(id);
    return hwnd != nullptr && SendMessageW(hwnd, CB_GETCURSEL, 0, 0) >= 0;
}

std::string get_edit_utf8(int id) {
    HWND hwnd = find_control(id);
    if (hwnd == nullptr) return {};
    return wide_to_utf8(get_window_text(hwnd));
}

void set_check(int id, bool checked) {
    HWND hwnd = find_control(id);
    if (hwnd != nullptr) SendMessageW(hwnd, BM_SETCHECK, checked ? BST_CHECKED : BST_UNCHECKED, 0);
}

bool get_check(int id) {
    HWND hwnd = find_control(id);
    return hwnd != nullptr && SendMessageW(hwnd, BM_GETCHECK, 0, 0) == BST_CHECKED;
}

void set_combo_text(int id, const wchar_t* value) {
    HWND hwnd = find_control(id);
    if (hwnd != nullptr) {
        const LRESULT index = SendMessageW(hwnd, CB_FINDSTRINGEXACT, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(value));
        if (index >= 0) SendMessageW(hwnd, CB_SETCURSEL, static_cast<WPARAM>(index), 0);
    }
}

void set_combo_items(int id, const std::vector<const wchar_t*>& items, const wchar_t* preferred) {
    HWND hwnd = find_control(id);
    if (hwnd == nullptr) return;

    SendMessageW(hwnd, CB_RESETCONTENT, 0, 0);
    for (const wchar_t* item : items) {
        SendMessageW(hwnd, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(item));
    }
    set_combo_text(id, preferred);
    if (!combo_has_selection(id) && !items.empty()) {
        SendMessageW(hwnd, CB_SETCURSEL, 0, 0);
    }
}

std::wstring get_combo_text(int id) {
    HWND hwnd = find_control(id);
    if (hwnd == nullptr) return {};
    const LRESULT index = SendMessageW(hwnd, CB_GETCURSEL, 0, 0);
    if (index < 0) return {};
    const LRESULT length = SendMessageW(hwnd, CB_GETLBTEXTLEN, static_cast<WPARAM>(index), 0);
    std::vector<wchar_t> buffer(static_cast<size_t>(length) + 1);
    SendMessageW(hwnd, CB_GETLBTEXT, static_cast<WPARAM>(index), reinterpret_cast<LPARAM>(buffer.data()));
    return std::wstring(buffer.data(), static_cast<size_t>(length));
}

bool read_double(int id, const wchar_t* label, double& value) {
    const std::string text = trim_copy(get_edit_utf8(id));
    try {
        size_t used = 0;
        value = std::stod(text, &used);
        if (used != text.size()) throw std::invalid_argument("trailing characters");
        return true;
    } catch (...) {
        MessageBoxW(g_hwnd, (std::wstring(L"Invalid number for ") + label).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
}

bool read_int(int id, const wchar_t* label, int& value) {
    const std::string text = trim_copy(get_edit_utf8(id));
    try {
        size_t used = 0;
        value = std::stoi(text, &used);
        if (used != text.size()) throw std::invalid_argument("trailing characters");
        return true;
    } catch (...) {
        MessageBoxW(g_hwnd, (std::wstring(L"Invalid integer for ") + label).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
}

bool combo_is(int id, const wchar_t* value);

double unit_scale_to_meters(int unit_id) {
    const std::wstring unit = get_combo_text(unit_id);
    if (unit == L"cm") return 1.0e-2;
    if (unit == L"mm") return 1.0e-3;
    if (unit == L"um") return 1.0e-6;
    return 1.0;
}

double unit_scale_to_seconds(int unit_id) {
    const std::wstring unit = get_combo_text(unit_id);
    if (unit == L"ms") return 1.0e-3;
    return 1.0;
}

double unit_scale_to_pascals(int unit_id) {
    const std::wstring unit = get_combo_text(unit_id);
    if (unit == L"kPa") return 1.0e3;
    if (unit == L"MPa") return 1.0e6;
    return 1.0;
}

bool read_length_meters(int id, int unit_id, const wchar_t* label, double& meters) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    meters = value * unit_scale_to_meters(unit_id);
    return true;
}

bool read_time_seconds(int id, int unit_id, const wchar_t* label, double& seconds) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    seconds = value * unit_scale_to_seconds(unit_id);
    return true;
}

bool read_pressure_pascals(int id, int unit_id, const wchar_t* label, double& pascals) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    pascals = value * unit_scale_to_pascals(unit_id);
    return true;
}

bool read_angle_degrees(int id, int unit_id, const wchar_t* label, double& degrees) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    degrees = combo_is(unit_id, L"rad") ? value * 180.0 / kPi : value;
    return true;
}

bool read_tilt_slope(int id, int unit_id, const wchar_t* label, double& slope) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    if (combo_is(unit_id, L"deg")) {
        slope = std::tan(value * kPi / 180.0);
    } else if (combo_is(unit_id, L"rad")) {
        slope = std::tan(value);
    } else {
        slope = value;
    }
    return true;
}

bool read_omega_rad_s(int id, int unit_id, const wchar_t* label, double& omega) {
    double value = 0.0;
    if (!read_double(id, label, value)) return false;
    omega = combo_is(unit_id, L"rpm") ? value * 2.0 * kPi / 60.0 : value;
    return true;
}

BCType bc_type_from_combo(int id) {
    const std::wstring value = get_combo_text(id);
    if (value == L"NEUMANN") return BCType::NEUMANN;
    if (value == L"INLET_OUTLET") return BCType::INLET_OUTLET;
    return BCType::DIRICHLET;
}

void set_bc_combo(int id, BCType type) {
    set_combo_text(id, utf8_to_wide(to_config_value(type)).c_str());
}

bool combo_is(int id, const wchar_t* value) {
    return get_combo_text(id) == value;
}

bool unit_is_native_length(int unit_id) {
    return combo_is(unit_id, L"m") || get_combo_text(unit_id).empty();
}

bool unit_is_native_time(int unit_id) {
    return combo_is(unit_id, L"s") || get_combo_text(unit_id).empty();
}

bool unit_is_native_pressure(int unit_id) {
    return combo_is(unit_id, L"Pa") || get_combo_text(unit_id).empty();
}

void set_native_unit_defaults() {
    if (!combo_has_selection(kComboUnitR)) set_combo_text(kComboUnitR, L"m");
    if (!combo_has_selection(kComboUnitC)) set_combo_text(kComboUnitC, L"m");
    if (!combo_has_selection(kComboUnitE)) set_combo_text(kComboUnitE, L"m");
    if (!combo_has_selection(kComboUnitL)) set_combo_text(kComboUnitL, L"m");
    if (!combo_has_selection(kComboUnitAttitude)) set_combo_text(kComboUnitAttitude, L"deg");
    if (!combo_has_selection(kComboUnitLoadAngle)) set_combo_text(kComboUnitLoadAngle, L"deg");
    if (!combo_has_selection(kComboUnitTiltX)) set_combo_text(kComboUnitTiltX, L"m/m");
    if (!combo_has_selection(kComboUnitTiltY)) set_combo_text(kComboUnitTiltY, L"m/m");
    if (!combo_has_selection(kComboUnitEndT)) set_combo_text(kComboUnitEndT, L"s");
    if (!combo_has_selection(kComboUnitDt)) set_combo_text(kComboUnitDt, L"s");
    if (!combo_has_selection(kComboUnitWriteInterval)) set_combo_text(kComboUnitWriteInterval, L"s");
    if (!combo_has_selection(kComboUnitOmega)) set_combo_text(kComboUnitOmega, L"rad/s");
    if (!combo_has_selection(kComboUnitInletTheta)) set_combo_text(kComboUnitInletTheta, L"deg");
    if (!combo_has_selection(kComboUnitInletZ)) set_combo_text(kComboUnitInletZ, L"m");
    if (!combo_has_selection(kComboUnitInletSize)) set_combo_text(kComboUnitInletSize, L"deg");
    if (!combo_has_selection(kComboUnitInletPressure)) set_combo_text(kComboUnitInletPressure, L"Pa");
    if (!combo_has_selection(kComboUnitNTheta)) set_combo_text(kComboUnitNTheta, L"cells");
    if (!combo_has_selection(kComboUnitNZ)) set_combo_text(kComboUnitNZ, L"cells");
    if (!combo_has_selection(kComboUnitMu)) set_combo_text(kComboUnitMu, L"Pa s");
    if (!combo_has_selection(kComboUnitRho)) set_combo_text(kComboUnitRho, L"kg/m3");
    if (!combo_has_selection(kComboUnitPCav)) set_combo_text(kComboUnitPCav, L"Pa");
    if (!combo_has_selection(kComboUnitBulkModulus)) set_combo_text(kComboUnitBulkModulus, L"Pa");
    if (!combo_has_selection(kComboUnitBcSouthVal)) set_combo_text(kComboUnitBcSouthVal, L"Pa");
    if (!combo_has_selection(kComboUnitBcSouthTheta)) set_combo_text(kComboUnitBcSouthTheta, L"-");
    if (!combo_has_selection(kComboUnitBcNorthVal)) set_combo_text(kComboUnitBcNorthVal, L"Pa");
    if (!combo_has_selection(kComboUnitBcNorthTheta)) set_combo_text(kComboUnitBcNorthTheta, L"-");
    if (!combo_has_selection(kComboUnitMaxOuter)) set_combo_text(kComboUnitMaxOuter, L"iter");
    if (!combo_has_selection(kComboUnitOuterTol)) set_combo_text(kComboUnitOuterTol, L"-");
    if (!combo_has_selection(kComboUnitThetaMin)) set_combo_text(kComboUnitThetaMin, L"-");
}

void set_length_from_config_text(int edit_id, int unit_id, const char* key, double meters) {
    if (unit_is_native_length(unit_id)) {
        set_edit_from_config_text(edit_id, key, meters);
    } else {
        set_edit_text(edit_id, compact_number(meters / unit_scale_to_meters(unit_id)));
    }
}

void set_time_from_config_text(int edit_id, int unit_id, const char* key, double seconds) {
    if (unit_is_native_time(unit_id)) {
        set_edit_from_config_text(edit_id, key, seconds);
    } else {
        set_edit_text(edit_id, compact_number(seconds / unit_scale_to_seconds(unit_id)));
    }
}

void set_pressure_from_config_text(int edit_id, int unit_id, const char* key, double pascals) {
    if (unit_is_native_pressure(unit_id)) {
        set_edit_from_config_text(edit_id, key, pascals);
    } else {
        set_edit_text(edit_id, compact_number(pascals / unit_scale_to_pascals(unit_id)));
    }
}

std::string native_or_converted_text(int edit_id, int unit_id, double native_value, bool native_unit) {
    return native_unit ? trim_copy(get_edit_utf8(edit_id)) : compact_number(native_value);
}

CavitationModel cavitation_from_combo() {
    const std::wstring value = get_combo_text(kComboCavitation);
    if (value == L"GUMBEL") return CavitationModel::GUMBEL;
    if (value == L"FULL_SOMMERFELD") return CavitationModel::FULL_SOMMERFELD;
    return CavitationModel::ELROD_ADAMS;
}

std::string inlets_to_text(const SimulationConfig& cfg) {
    std::ostringstream out;
    for (const auto& inlet : cfg.inlets) {
        if (inlet.type == InletConfig::Type::CIRCULAR) {
            out << "inlet_circular = " << compact_number(inlet.theta) << " " << compact_number(inlet.z) << " "
                << compact_number(inlet.size) << " " << compact_number(inlet.p_supply) << "\r\n";
        } else {
            out << "inlet_groove = " << compact_number(inlet.theta) << " "
                << compact_number(inlet.size) << " " << compact_number(inlet.p_supply) << "\r\n";
        }
    }
    return out.str();
}

struct InletTextFields {
    bool valid = false;
    InletConfig::Type type = InletConfig::Type::GROOVE;
    std::string theta;
    std::string z;
    std::string size;
    std::string p_supply;
};

InletTextFields first_inlet_text_fields(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        const size_t comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);
        line = trim_copy(line);
        if (line.empty()) continue;

        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = normalise_key(line.substr(0, eq));
        std::stringstream values(line.substr(eq + 1));

        InletTextFields fields;
        if (key == "inlet_circular") {
            fields.type = InletConfig::Type::CIRCULAR;
            if (values >> fields.theta >> fields.z >> fields.size >> fields.p_supply) fields.valid = true;
        } else if (key == "inlet_groove") {
            fields.type = InletConfig::Type::GROOVE;
            if (values >> fields.theta >> fields.size >> fields.p_supply) fields.valid = true;
        }
        if (fields.valid) return fields;
    }
    return {};
}

void update_inlet_editor_visibility();
void layout_parameter_panel();

void populate_inlet_editor(const SimulationConfig& cfg) {
    InletConfig inlet;
    if (!cfg.inlets.empty()) {
        inlet = cfg.inlets.front();
    } else {
        inlet.type = InletConfig::Type::GROOVE;
        inlet.theta = 90.0;
        inlet.z = 0.5 * cfg.L;
        inlet.size = 10.0;
        inlet.p_supply = 2e5;
    }

    set_combo_text(kComboInletType, inlet.type == InletConfig::Type::CIRCULAR ? L"CIRCULAR" : L"GROOVE");
    update_inlet_editor_visibility();
    const InletTextFields raw = first_inlet_text_fields(g_config_inlets_text);
    const bool use_raw = raw.valid && raw.type == inlet.type;
    if (combo_is(kComboUnitInletTheta, L"rad")) {
        set_edit_text(kEditInletTheta, compact_number(inlet.theta * kPi / 180.0));
    } else {
        set_edit_text(kEditInletTheta, use_raw ? raw.theta : compact_number(inlet.theta));
    }
    set_edit_text(kEditInletZ, (use_raw && unit_is_native_length(kComboUnitInletZ)) ? raw.z : compact_number(inlet.z / unit_scale_to_meters(kComboUnitInletZ)));
    if (inlet.type == InletConfig::Type::CIRCULAR) {
        set_edit_text(kEditInletSize, (use_raw && unit_is_native_length(kComboUnitInletSize)) ? raw.size : compact_number(inlet.size / unit_scale_to_meters(kComboUnitInletSize)));
    } else if (combo_is(kComboUnitInletSize, L"rad")) {
        set_edit_text(kEditInletSize, compact_number(inlet.size * kPi / 180.0));
    } else {
        set_edit_text(kEditInletSize, use_raw ? raw.size : compact_number(inlet.size));
    }
    set_edit_text(
        kEditInletPressure,
        (use_raw && unit_is_native_pressure(kComboUnitInletPressure))
            ? raw.p_supply
            : compact_number(inlet.p_supply / unit_scale_to_pascals(kComboUnitInletPressure)));
}

void update_inlet_editor_visibility() {
    const bool circular = combo_is(kComboInletType, L"CIRCULAR");
    ShowWindow(g_inlet_z_label, circular ? SW_SHOW : SW_HIDE);
    ShowWindow(find_control(kEditInletZ), circular ? SW_SHOW : SW_HIDE);
    ShowWindow(find_control(kComboUnitInletZ), circular ? SW_SHOW : SW_HIDE);
    set_combo_items(
        kComboUnitInletSize,
        circular ? std::vector<const wchar_t*>{L"m", L"cm", L"mm", L"um"} : std::vector<const wchar_t*>{L"deg", L"rad"},
        circular ? L"m" : L"deg");
    if (g_inlet_size_label != nullptr) {
        SetWindowTextW(g_inlet_size_label, circular ? L"radius" : L"width");
    }
    layout_parameter_panel();
}

bool read_inlet_editor(InletConfig& inlet) {
    inlet.type = combo_is(kComboInletType, L"CIRCULAR") ? InletConfig::Type::CIRCULAR : InletConfig::Type::GROOVE;
    if (!read_angle_degrees(kEditInletTheta, kComboUnitInletTheta, L"inlet theta", inlet.theta)) return false;
    if (inlet.type == InletConfig::Type::CIRCULAR && !read_length_meters(kEditInletZ, kComboUnitInletZ, L"inlet z", inlet.z)) return false;
    if (inlet.type == InletConfig::Type::CIRCULAR) {
        if (!read_length_meters(kEditInletSize, kComboUnitInletSize, L"inlet radius", inlet.size)) return false;
    } else if (!read_angle_degrees(kEditInletSize, kComboUnitInletSize, L"inlet width", inlet.size)) {
        return false;
    }
    if (!read_pressure_pascals(kEditInletPressure, kComboUnitInletPressure, L"inlet p_supply", inlet.p_supply)) return false;
    if (inlet.type == InletConfig::Type::GROOVE) inlet.z = 0.0;
    return true;
}

bool parse_inlets(const std::string& text, std::vector<InletConfig>& inlets) {
    inlets.clear();
    std::istringstream input(text);
    std::string line;
    int line_number = 0;

    while (std::getline(input, line)) {
        ++line_number;
        const size_t comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);
        line = trim_copy(line);
        if (line.empty()) continue;

        const size_t eq = line.find('=');
        if (eq == std::string::npos) {
            MessageBoxW(
                g_hwnd,
                (std::wstring(L"Inlet line ") + std::to_wstring(line_number) + L" is missing '='.").c_str(),
                L"Pancake",
                MB_ICONERROR | MB_OK);
            return false;
        }

        const std::string key = normalise_key(line.substr(0, eq));
        std::stringstream values(line.substr(eq + 1));
        InletConfig inlet;

        if (key == "inlet_circular") {
            inlet.type = InletConfig::Type::CIRCULAR;
            if (!(values >> inlet.theta >> inlet.z >> inlet.size >> inlet.p_supply)) {
                MessageBoxW(g_hwnd, L"Expected inlet_circular = theta z radius p_supply.", L"Pancake", MB_ICONERROR | MB_OK);
                return false;
            }
            inlets.push_back(inlet);
        } else if (key == "inlet_groove") {
            inlet.type = InletConfig::Type::GROOVE;
            if (!(values >> inlet.theta >> inlet.size >> inlet.p_supply)) {
                MessageBoxW(g_hwnd, L"Expected inlet_groove = theta width p_supply.", L"Pancake", MB_ICONERROR | MB_OK);
                return false;
            }
            inlets.push_back(inlet);
        } else {
            MessageBoxW(
                g_hwnd,
                (std::wstring(L"Unknown inlet key on line ") + std::to_wstring(line_number) + L".").c_str(),
                L"Pancake",
                MB_ICONERROR | MB_OK);
            return false;
        }
    }

    return true;
}

void populate_preview_fields() {
    HWND combo = find_control(kComboPreviewField);
    if (combo == nullptr) return;

    SendMessageW(combo, CB_RESETCONTENT, 0, 0);
    for (const auto& field : kOutputFields) {
        if (std::string(field.key) == "velocity") continue;
        SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(field.label));
    }
    SendMessageW(combo, CB_SETCURSEL, 0, 0);
}

std::string form_config_text(const SimulationConfig& cfg);

void populate_form_from_config(const SimulationConfig& cfg) {
    g_loading = true;

    set_native_unit_defaults();

    set_length_from_config_text(kEditR, kComboUnitR, "R", cfg.R);
    set_length_from_config_text(kEditC, kComboUnitC, "c", cfg.c);
    set_length_from_config_text(kEditE, kComboUnitE, "e", cfg.e);
    set_length_from_config_text(kEditL, kComboUnitL, "L", cfg.L);
    if (combo_is(kComboUnitAttitude, L"rad")) {
        set_edit_text(kEditAttitude, compact_number(cfg.attitude_angle_deg * kPi / 180.0));
    } else {
        set_edit_from_config_text(kEditAttitude, "attitude_angle_deg", cfg.attitude_angle_deg);
    }
    if (combo_is(kComboUnitLoadAngle, L"rad")) {
        set_edit_text(kEditLoadAngle, compact_number(cfg.load_angle_deg * kPi / 180.0));
    } else {
        set_edit_from_config_text(kEditLoadAngle, "load_angle_deg", cfg.load_angle_deg);
    }
    if (combo_is(kComboUnitTiltX, L"deg")) {
        set_edit_text(kEditTiltX, compact_number(std::atan(cfg.tilt_slope_x) * 180.0 / kPi));
    } else if (combo_is(kComboUnitTiltX, L"rad")) {
        set_edit_text(kEditTiltX, compact_number(std::atan(cfg.tilt_slope_x)));
    } else {
        set_edit_from_config_text(kEditTiltX, "tilt_slope_x", cfg.tilt_slope_x);
    }
    if (combo_is(kComboUnitTiltY, L"deg")) {
        set_edit_text(kEditTiltY, compact_number(std::atan(cfg.tilt_slope_y) * 180.0 / kPi));
    } else if (combo_is(kComboUnitTiltY, L"rad")) {
        set_edit_text(kEditTiltY, compact_number(std::atan(cfg.tilt_slope_y)));
    } else {
        set_edit_from_config_text(kEditTiltY, "tilt_slope_y", cfg.tilt_slope_y);
    }
    set_edit_from_config_text(kEditNTheta, "n_theta_global", cfg.n_theta_global);
    set_edit_from_config_text(kEditNZ, "n_z_global", cfg.n_z_global);
    set_time_from_config_text(kEditEndT, kComboUnitEndT, "end_t", cfg.end_t);
    set_time_from_config_text(kEditDt, kComboUnitDt, "dt", cfg.dt);
    set_time_from_config_text(kEditWriteInterval, kComboUnitWriteInterval, "write_interval", cfg.write_interval);
    if (combo_is(kComboUnitOmega, L"rpm")) {
        set_edit_text(kEditOmega, compact_number(cfg.omega * 60.0 / (2.0 * kPi)));
    } else {
        set_edit_from_config_text(kEditOmega, "omega", cfg.omega);
    }
    set_edit_from_config_text(kEditMu, "mu", cfg.mu);
    set_edit_from_config_text(kEditRho, "rho", cfg.rho);
    set_pressure_from_config_text(kEditPCav, kComboUnitPCav, "p_cav", cfg.p_cav);
    set_pressure_from_config_text(kEditBulkModulus, kComboUnitBulkModulus, "bulk_modulus", cfg.bulk_modulus);
    set_combo_text(kComboCavitation, utf8_to_wide(to_config_value(cfg.cavitation_model)).c_str());
    set_edit_from_config_text(kEditMaxOuter, "max_outer_iters", cfg.max_outer_iters);
    set_edit_from_config_text(kEditOuterTol, "outer_tol", cfg.outer_tol);
    set_edit_from_config_text(kEditThetaMin, "theta_min", cfg.theta_min);
    set_check(kCheckLogOuter, cfg.log_outer_iters);
    set_bc_combo(kComboBcSouth, cfg.bc_z_south_type);
    set_pressure_from_config_text(kEditBcSouthVal, kComboUnitBcSouthVal, "bc_z_south_val", cfg.bc_z_south_val);
    set_edit_from_config_text(kEditBcSouthTheta, "bc_z_south_theta", cfg.bc_z_south_theta);
    set_bc_combo(kComboBcNorth, cfg.bc_z_north_type);
    set_pressure_from_config_text(kEditBcNorthVal, kComboUnitBcNorthVal, "bc_z_north_val", cfg.bc_z_north_val);
    set_edit_from_config_text(kEditBcNorthTheta, "bc_z_north_theta", cfg.bc_z_north_theta);
    populate_inlet_editor(cfg);
    update_inlet_editor_visibility();

    set_edit_text(kEditOutputDir, cfg.output_dir);
    set_edit_text(kEditFilenamePrefix, cfg.filename_prefix);
    set_check(kCheckOutput3d, cfg.output_write_3d);
    set_check(kCheckOutputFlat, cfg.output_write_flat);
    for (const auto& field : kOutputFields) {
        set_check(field.id, cfg.output_field_enabled(field.key));
    }

    const std::string current_solver_path = trim_copy(get_edit_utf8(kEditSolverPath));
    const std::string current_ranks = trim_copy(get_edit_utf8(kEditMpiRanks));
    if (current_solver_path.empty()) {
        const fs::path solver_path = g_executable_dir / L"pancake.exe";
        set_edit_text(kEditSolverPath, wide_to_utf8(solver_path.wstring()));
    }
    if (current_ranks.empty()) {
        set_edit_text(kEditMpiRanks, 1);
    }
    set_edit_text(kEditRawConfig, form_config_text(cfg));

    g_loading = false;
}

bool apply_form_to_config(SimulationConfig& cfg) {
    if (!read_length_meters(kEditR, kComboUnitR, L"R", cfg.R)) return false;
    if (!read_length_meters(kEditC, kComboUnitC, L"c", cfg.c)) return false;
    if (!read_length_meters(kEditE, kComboUnitE, L"e", cfg.e)) return false;
    if (!read_length_meters(kEditL, kComboUnitL, L"L", cfg.L)) return false;
    if (!read_angle_degrees(kEditAttitude, kComboUnitAttitude, L"attitude_angle", cfg.attitude_angle_deg)) return false;
    if (!read_angle_degrees(kEditLoadAngle, kComboUnitLoadAngle, L"load_angle", cfg.load_angle_deg)) return false;
    if (!read_tilt_slope(kEditTiltX, kComboUnitTiltX, L"tilt_slope_x", cfg.tilt_slope_x)) return false;
    if (!read_tilt_slope(kEditTiltY, kComboUnitTiltY, L"tilt_slope_y", cfg.tilt_slope_y)) return false;
    if (!read_int(kEditNTheta, L"n_theta_global", cfg.n_theta_global)) return false;
    if (!read_int(kEditNZ, L"n_z_global", cfg.n_z_global)) return false;
    if (!read_time_seconds(kEditEndT, kComboUnitEndT, L"end_t", cfg.end_t)) return false;
    if (!read_time_seconds(kEditDt, kComboUnitDt, L"dt", cfg.dt)) return false;
    if (!read_time_seconds(kEditWriteInterval, kComboUnitWriteInterval, L"write_interval", cfg.write_interval)) return false;
    if (!read_omega_rad_s(kEditOmega, kComboUnitOmega, L"omega", cfg.omega)) return false;
    if (!read_double(kEditMu, L"mu", cfg.mu)) return false;
    if (!read_double(kEditRho, L"rho", cfg.rho)) return false;
    if (!read_pressure_pascals(kEditPCav, kComboUnitPCav, L"p_cav", cfg.p_cav)) return false;
    if (!read_pressure_pascals(kEditBulkModulus, kComboUnitBulkModulus, L"bulk_modulus", cfg.bulk_modulus)) return false;
    cfg.cavitation_model = cavitation_from_combo();
    if (!read_int(kEditMaxOuter, L"max_outer_iters", cfg.max_outer_iters)) return false;
    if (!read_double(kEditOuterTol, L"outer_tol", cfg.outer_tol)) return false;
    if (!read_double(kEditThetaMin, L"theta_min", cfg.theta_min)) return false;
    cfg.log_outer_iters = get_check(kCheckLogOuter);

    cfg.bc_z_south_type = bc_type_from_combo(kComboBcSouth);
    if (!read_pressure_pascals(kEditBcSouthVal, kComboUnitBcSouthVal, L"bc_z_south_val", cfg.bc_z_south_val)) return false;
    if (!read_double(kEditBcSouthTheta, L"bc_z_south_theta", cfg.bc_z_south_theta)) return false;
    cfg.bc_z_north_type = bc_type_from_combo(kComboBcNorth);
    if (!read_pressure_pascals(kEditBcNorthVal, kComboUnitBcNorthVal, L"bc_z_north_val", cfg.bc_z_north_val)) return false;
    if (!read_double(kEditBcNorthTheta, L"bc_z_north_theta", cfg.bc_z_north_theta)) return false;

    InletConfig inlet;
    if (!read_inlet_editor(inlet)) return false;
    cfg.inlets = {inlet};

    cfg.output_dir = trim_copy(get_edit_utf8(kEditOutputDir));
    cfg.filename_prefix = trim_copy(get_edit_utf8(kEditFilenamePrefix));
    cfg.output_write_3d = get_check(kCheckOutput3d);
    cfg.output_write_flat = get_check(kCheckOutputFlat);
    cfg.output_fields.clear();
    for (const auto& field : kOutputFields) {
        if (get_check(field.id)) cfg.output_fields.push_back(field.key);
    }

    if (cfg.output_dir.empty()) {
        MessageBoxW(g_hwnd, L"output_dir cannot be empty.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    if (cfg.filename_prefix.empty()) {
        MessageBoxW(g_hwnd, L"filename_prefix cannot be empty.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    if (cfg.n_theta_global <= 0 || cfg.n_z_global <= 0) {
        MessageBoxW(g_hwnd, L"Grid sizes must be positive.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    if (cfg.dt <= 0.0 || cfg.write_interval <= 0.0) {
        MessageBoxW(g_hwnd, L"dt and write_interval must be positive.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }

    return true;
}

std::string form_config_text(const SimulationConfig& cfg) {
    std::ostringstream out;
    out << "# Geometry\n"
        << "R = " << native_or_converted_text(kEditR, kComboUnitR, cfg.R, unit_is_native_length(kComboUnitR)) << "\n"
        << "c = " << native_or_converted_text(kEditC, kComboUnitC, cfg.c, unit_is_native_length(kComboUnitC)) << "\n"
        << "e = " << native_or_converted_text(kEditE, kComboUnitE, cfg.e, unit_is_native_length(kComboUnitE)) << "\n"
        << "L = " << native_or_converted_text(kEditL, kComboUnitL, cfg.L, unit_is_native_length(kComboUnitL)) << "\n"
        << "attitude_angle_deg = " << (combo_is(kComboUnitAttitude, L"deg") ? trim_copy(get_edit_utf8(kEditAttitude)) : compact_number(cfg.attitude_angle_deg)) << "\n"
        << "load_angle_deg = " << (combo_is(kComboUnitLoadAngle, L"deg") ? trim_copy(get_edit_utf8(kEditLoadAngle)) : compact_number(cfg.load_angle_deg)) << "\n"
        << "tilt_slope_x = " << (combo_is(kComboUnitTiltX, L"m/m") ? trim_copy(get_edit_utf8(kEditTiltX)) : compact_number(cfg.tilt_slope_x)) << "\n"
        << "tilt_slope_y = " << (combo_is(kComboUnitTiltY, L"m/m") ? trim_copy(get_edit_utf8(kEditTiltY)) : compact_number(cfg.tilt_slope_y)) << "\n\n";

    out << "# Grid\n"
        << "n_theta_global = " << trim_copy(get_edit_utf8(kEditNTheta)) << "\n"
        << "n_z_global = " << trim_copy(get_edit_utf8(kEditNZ)) << "\n\n";

    out << "# Time\n"
        << "end_t = " << native_or_converted_text(kEditEndT, kComboUnitEndT, cfg.end_t, unit_is_native_time(kComboUnitEndT)) << "\n"
        << "dt = " << native_or_converted_text(kEditDt, kComboUnitDt, cfg.dt, unit_is_native_time(kComboUnitDt)) << "\n"
        << "write_interval = " << native_or_converted_text(kEditWriteInterval, kComboUnitWriteInterval, cfg.write_interval, unit_is_native_time(kComboUnitWriteInterval)) << "\n\n";

    out << "# Physics\n"
        << "omega = " << (combo_is(kComboUnitOmega, L"rad/s") ? trim_copy(get_edit_utf8(kEditOmega)) : compact_number(cfg.omega)) << "\n"
        << "mu = " << trim_copy(get_edit_utf8(kEditMu)) << "\n"
        << "rho = " << trim_copy(get_edit_utf8(kEditRho)) << "\n"
        << "p_cav = " << native_or_converted_text(kEditPCav, kComboUnitPCav, cfg.p_cav, unit_is_native_pressure(kComboUnitPCav)) << "\n"
        << "bulk_modulus = " << native_or_converted_text(kEditBulkModulus, kComboUnitBulkModulus, cfg.bulk_modulus, unit_is_native_pressure(kComboUnitBulkModulus)) << "\n\n";

    out << "# Axial Boundary Conditions (Outlets)\n"
        << "# Options: DIRICHLET, NEUMANN, INLET_OUTLET\n"
        << "bc_z_south_type = " << to_config_value(cfg.bc_z_south_type) << "\n"
        << "bc_z_south_val = " << native_or_converted_text(kEditBcSouthVal, kComboUnitBcSouthVal, cfg.bc_z_south_val, unit_is_native_pressure(kComboUnitBcSouthVal)) << "\n"
        << "bc_z_south_theta = " << trim_copy(get_edit_utf8(kEditBcSouthTheta)) << "\n"
        << "bc_z_north_type = " << to_config_value(cfg.bc_z_north_type) << "\n"
        << "bc_z_north_val = " << native_or_converted_text(kEditBcNorthVal, kComboUnitBcNorthVal, cfg.bc_z_north_val, unit_is_native_pressure(kComboUnitBcNorthVal)) << "\n"
        << "bc_z_north_theta = " << trim_copy(get_edit_utf8(kEditBcNorthTheta)) << "\n\n";

    out << "# Cavitation\n"
        << "# Options: FULL_SOMMERFELD, GUMBEL, ELROD_ADAMS\n"
        << "cavitation_model = " << to_config_value(cfg.cavitation_model) << "\n"
        << "max_outer_iters = " << trim_copy(get_edit_utf8(kEditMaxOuter)) << "\n"
        << "outer_tol = " << trim_copy(get_edit_utf8(kEditOuterTol)) << "\n"
        << "theta_min = " << trim_copy(get_edit_utf8(kEditThetaMin)) << "\n"
        << "log_outer_iters = " << (cfg.log_outer_iters ? "true" : "false") << "\n\n";

    out << "# Inlets\n"
        << "# inlet_circular = theta(deg) z(m) radius(m) p_supply(Pa)\n"
        << "# inlet_groove   = theta(deg) width(deg) p_supply(Pa)\n";
    const InletConfig inlet = cfg.inlets.empty() ? InletConfig{} : cfg.inlets.front();
    const std::string inlet_theta =
        combo_is(kComboUnitInletTheta, L"deg") ? trim_copy(get_edit_utf8(kEditInletTheta)) : compact_number(inlet.theta);
    const std::string inlet_z =
        unit_is_native_length(kComboUnitInletZ) ? trim_copy(get_edit_utf8(kEditInletZ)) : compact_number(inlet.z);
    const std::string inlet_size = (inlet.type == InletConfig::Type::CIRCULAR)
        ? native_or_converted_text(kEditInletSize, kComboUnitInletSize, inlet.size, unit_is_native_length(kComboUnitInletSize))
        : (combo_is(kComboUnitInletSize, L"deg") ? trim_copy(get_edit_utf8(kEditInletSize)) : compact_number(inlet.size));
    const std::string inlet_pressure =
        native_or_converted_text(kEditInletPressure, kComboUnitInletPressure, inlet.p_supply, unit_is_native_pressure(kComboUnitInletPressure));
    if (!cfg.inlets.empty() && cfg.inlets.front().type == InletConfig::Type::CIRCULAR) {
        out << "inlet_circular = "
            << inlet_theta << " "
            << inlet_z << " "
            << inlet_size << " "
            << inlet_pressure << "\n"
            << "# inlet_groove = "
            << inlet_theta << " "
            << compact_number(inlet.size) << " "
            << inlet_pressure << "\n";
    } else {
        out << "inlet_groove = "
            << inlet_theta << " "
            << inlet_size << " "
            << inlet_pressure << "\n"
            << "# inlet_circular = "
            << inlet_theta << " "
            << compact_number(0.5 * cfg.L) << " "
            << compact_number(cfg.c) << " "
            << inlet_pressure << "\n";
    }
    out << "\n";

    out << "# Output\n"
        << "# output_fields accepts comma-separated field names.\n"
        << "# Options: pressure, film_content, h, rho, inlet_indicator, velocity,\n"
        << "#          load_x, load_y, load_z, friction_torque\n"
        << "output_dir = " << trim_copy(get_edit_utf8(kEditOutputDir)) << "\n"
        << "filename_prefix = " << trim_copy(get_edit_utf8(kEditFilenamePrefix)) << "\n"
        << "output_write_3d = " << (cfg.output_write_3d ? "true" : "false") << "\n"
        << "output_write_flat = " << (cfg.output_write_flat ? "true" : "false") << "\n"
        << "output_fields = " << cfg.output_fields_text() << "\n";

    return out.str();
}

bool write_text_file(const fs::path& path, const std::string& text) {
    std::ofstream file(path, std::ios::trunc);
    if (!file.is_open()) {
        MessageBoxW(g_hwnd, (std::wstring(L"Could not open config file for writing: ") + path.wstring()).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    file << text;
    return true;
}

bool save_config() {
    SimulationConfig cfg = g_config;
    std::string config_text;
    if (g_active_tab == kRawTab) {
        config_text = get_edit_utf8(kEditRawConfig);
        cfg = SimulationConfig{};
        try {
            cfg.load_from_text(config_text);
        } catch (const std::exception& ex) {
            MessageBoxW(g_hwnd, utf8_to_wide(ex.what()).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
            return false;
        }
    } else if (!apply_form_to_config(cfg)) {
        return false;
    } else {
        config_text = form_config_text(cfg);
    }

    if (!write_text_file(g_config_path, config_text)) return false;

    g_config = cfg;
    g_config_value_text = parse_config_value_text(config_text);
    g_config_inlets_text = parse_config_inlets_text(config_text);
    populate_form_from_config(g_config);
    g_dirty = false;
    set_status(std::wstring(L"Saved ") + g_config_path.wstring());
    return true;
}

bool confirm_discard_changes() {
    if (!g_dirty) return true;
    const int result = MessageBoxW(
        g_hwnd,
        L"Discard the unsaved configuration changes?",
        L"Pancake",
        MB_ICONWARNING | MB_YESNO | MB_DEFBUTTON2);
    return result == IDYES;
}

void load_config() {
    g_config = SimulationConfig{};
    g_config_value_text.clear();
    g_config_inlets_text.clear();
    if (fs::exists(g_config_path)) {
        const std::string text = read_file_text(g_config_path);
        g_config.load_from_text(text);
        g_config_value_text = parse_config_value_text(text);
        g_config_inlets_text = parse_config_inlets_text(text);
        set_status(std::wstring(L"Loaded ") + g_config_path.wstring());
    } else {
        set_status(std::wstring(L"Using defaults. Save creates ") + g_config_path.wstring());
    }

    populate_form_from_config(g_config);
    g_dirty = false;
}

void reset_defaults() {
    if (!confirm_discard_changes()) return;
    g_config = SimulationConfig{};
    g_config_value_text.clear();
    g_config_inlets_text.clear();
    populate_form_from_config(g_config);
    g_dirty = true;
    set_status(L"Loaded defaults. Save to write config.txt.");
}

fs::path output_directory_for_config(const SimulationConfig& cfg) {
    fs::path output_path = utf8_to_wide(cfg.output_dir);
    if (output_path.is_relative()) output_path = g_executable_dir / output_path;
    return output_path;
}

bool try_parse_ansi_codes(const std::wstring& text, size_t offset, size_t& next, std::vector<int>& codes) {
    size_t start = std::wstring::npos;
    if (text[offset] == L'\x1b' && offset + 1 < text.size() && text[offset + 1] == L'[') {
        start = offset + 2;
    } else if (text[offset] == L'[' && offset + 2 < text.size() && std::iswdigit(text[offset + 1])) {
        start = offset + 1;
    } else {
        return false;
    }

    size_t end = start;
    while (end < text.size() && end - start <= 16) {
        const wchar_t ch = text[end];
        if (ch == L'm') break;
        if (!std::iswdigit(ch) && ch != L';') return false;
        ++end;
    }
    if (end >= text.size() || text[end] != L'm' || end == start) return false;

    std::wstring code_text = text.substr(start, end - start);
    codes.clear();
    size_t code_start = 0;
    while (code_start <= code_text.size()) {
        const size_t code_end = code_text.find(L';', code_start);
        const std::wstring token = code_text.substr(code_start, code_end == std::wstring::npos ? std::wstring::npos : code_end - code_start);
        if (!token.empty()) {
            try {
                codes.push_back(std::stoi(token));
            } catch (...) {
                return false;
            }
        }
        if (code_end == std::wstring::npos) break;
        code_start = code_end + 1;
    }

    next = end + 1;
    return !codes.empty();
}

void apply_ansi_codes(const std::vector<int>& codes, COLORREF& color) {
    for (int code : codes) {
        switch (code) {
            case 0:
            case 39:
                color = RGB(38, 46, 55);
                break;
            case 31:
            case 91:
                color = RGB(183, 61, 69);
                break;
            case 32:
            case 92:
                color = RGB(25, 128, 74);
                break;
            case 33:
            case 93:
                color = RGB(153, 109, 21);
                break;
            case 34:
            case 94:
                color = RGB(40, 110, 196);
                break;
            case 35:
            case 95:
                color = RGB(150, 70, 168);
                break;
            case 36:
            case 96:
                color = RGB(20, 119, 145);
                break;
            case 90:
                color = RGB(99, 110, 123);
                break;
            default:
                break;
        }
    }
}

void append_log_segment(HWND log, const std::wstring& segment, COLORREF color) {
    if (segment.empty()) return;
    const std::wstring normalized = normalize_newlines_for_edit(segment);
    const int length = GetWindowTextLengthW(log);
    SendMessageW(log, EM_SETSEL, static_cast<WPARAM>(length), static_cast<LPARAM>(length));
    if (g_richedit_available) {
        CHARFORMAT2W format{};
        format.cbSize = sizeof(format);
        format.dwMask = CFM_COLOR;
        format.crTextColor = color;
        SendMessageW(log, EM_SETCHARFORMAT, SCF_SELECTION, reinterpret_cast<LPARAM>(&format));
    }
    SendMessageW(log, EM_REPLACESEL, FALSE, reinterpret_cast<LPARAM>(normalized.c_str()));
}

void append_log(const std::wstring& text) {
    HWND log = find_control(kEditLog);
    if (log == nullptr) return;

    COLORREF color = RGB(38, 46, 55);
    size_t segment_start = 0;
    size_t i = 0;
    while (i < text.size()) {
        size_t next = 0;
        std::vector<int> codes;
        if (try_parse_ansi_codes(text, i, next, codes)) {
            append_log_segment(log, text.substr(segment_start, i - segment_start), color);
            apply_ansi_codes(codes, color);
            i = next;
            segment_start = i;
        } else {
            ++i;
        }
    }
    append_log_segment(log, text.substr(segment_start), color);
    SendMessageW(log, EM_SCROLLCARET, 0, 0);
}

void clear_log() {
    HWND log = find_control(kEditLog);
    if (log != nullptr) SetWindowTextW(log, L"");
}

std::wstring quote_arg(const std::wstring& arg) {
    std::wstring result = L"\"";
    for (wchar_t ch : arg) {
        if (ch == L'"') result += L'\\';
        result += ch;
    }
    result += L"\"";
    return result;
}

std::wstring get_environment_value(const wchar_t* name) {
    const DWORD length = GetEnvironmentVariableW(name, nullptr, 0);
    if (length == 0) return {};

    std::wstring value(length, L'\0');
    const DWORD written = GetEnvironmentVariableW(name, value.data(), length);
    if (written == 0) return {};
    value.resize(written);
    return value;
}

fs::path find_mpiexec() {
    std::vector<fs::path> candidates;

    const std::wstring msmpi_bin = get_environment_value(L"MSMPI_BIN");
    if (!msmpi_bin.empty()) candidates.push_back(fs::path(msmpi_bin) / L"mpiexec.exe");

    candidates.push_back(fs::path(L"C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"));
    candidates.push_back(fs::path(L"C:/Program Files (x86)/Microsoft SDKs/MPI/Bin/mpiexec.exe"));

    for (const auto& candidate : candidates) {
        if (fs::exists(candidate)) return candidate;
    }

    wchar_t path[MAX_PATH] = {};
    const DWORD length = SearchPathW(nullptr, L"mpiexec.exe", nullptr, MAX_PATH, path, nullptr);
    if (length > 0 && length < MAX_PATH) return fs::path(path);

    return {};
}

void close_solver_handles() {
    if (g_solver.stdout_read != nullptr) {
        CloseHandle(g_solver.stdout_read);
        g_solver.stdout_read = nullptr;
    }
    if (g_solver.process_info.hThread != nullptr) {
        CloseHandle(g_solver.process_info.hThread);
        g_solver.process_info.hThread = nullptr;
    }
    if (g_solver.process_info.hProcess != nullptr) {
        CloseHandle(g_solver.process_info.hProcess);
        g_solver.process_info.hProcess = nullptr;
    }
    g_solver.process_info.dwProcessId = 0;
    g_solver.running = false;
}

void set_run_state(bool running) {
    EnableWindow(g_run_button, running ? FALSE : TRUE);
    EnableWindow(g_stop_button, running ? TRUE : FALSE);
    if (g_run_panel_button != nullptr) EnableWindow(g_run_panel_button, running ? FALSE : TRUE);
    if (g_stop_panel_button != nullptr) EnableWindow(g_stop_panel_button, running ? TRUE : FALSE);
}

bool read_mpi_ranks(int& ranks) {
    if (!read_int(kEditMpiRanks, L"MPI ranks", ranks)) return false;
    if (ranks < 1) {
        MessageBoxW(g_hwnd, L"MPI ranks must be at least 1.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    return true;
}

fs::path project_root_from_gui() {
    const fs::path parent = g_executable_dir.parent_path();
    if (fs::exists(parent / L"CMakeLists.txt")) return parent;
    return fs::current_path();
}

bool make_solver_command(int ranks, std::wstring& command, fs::path& working_dir, std::wstring& missing_message) {
    const fs::path project_root = project_root_from_gui();
    const std::wstring requested_text = utf8_to_wide(trim_copy(get_edit_utf8(kEditSolverPath)));

    std::vector<fs::path> native_candidates;
    if (!requested_text.empty()) {
        fs::path requested = requested_text;
        if (requested.is_relative()) requested = g_executable_dir / requested;
        native_candidates.push_back(requested);
    }
    native_candidates.push_back(g_executable_dir / L"pancake.exe");
    native_candidates.push_back(project_root / L"build-windows-mingw" / L"pancake.exe");

    for (const fs::path& solver_path : native_candidates) {
        if (solver_path.empty() || !fs::exists(solver_path)) continue;
        if (ranks > 1) {
            const fs::path mpiexec_path = find_mpiexec();
            if (mpiexec_path.empty()) {
                missing_message = L"Microsoft MPI mpiexec.exe was not found. Install the Microsoft MPI Runtime, then run again.";
                return false;
            }
            command = quote_arg(mpiexec_path.wstring()) + L" -n " + std::to_wstring(ranks) + L" " +
                      quote_arg(solver_path.wstring()) + L" -c " + quote_arg(g_config_path.wstring());
        } else {
            command = quote_arg(solver_path.wstring()) + L" -c " + quote_arg(g_config_path.wstring());
        }
        working_dir = g_executable_dir;
        return true;
    }

    missing_message =
        L"No solver backend was found.\n\n"
        L"Build the native Windows workflow from MSYS2 MINGW64:\n\n"
        L"  cmake --fresh --preset windows-native-mingw\n"
        L"  cmake --build --preset windows-native-mingw-release\n\n"
        L"The GUI expects pancake.exe beside pancake_gui.exe or under build-windows-mingw.";
    return false;
}

void start_solver() {
    if (g_solver.running) return;
    if (!save_config()) return;

    int ranks = 1;
    if (!read_mpi_ranks(ranks)) return;

    std::wstring command;
    std::wstring missing_message;
    fs::path working_dir = g_executable_dir;
    if (!make_solver_command(ranks, command, working_dir, missing_message)) {
        MessageBoxW(g_hwnd, missing_message.c_str(), L"Pancake", MB_ICONINFORMATION | MB_OK);
        set_status(L"No solver backend found.");
        return;
    }

    SECURITY_ATTRIBUTES security{};
    security.nLength = sizeof(security);
    security.bInheritHandle = TRUE;

    HANDLE read_pipe = nullptr;
    HANDLE write_pipe = nullptr;
    if (!CreatePipe(&read_pipe, &write_pipe, &security, 0)) {
        MessageBoxW(g_hwnd, L"Could not create solver output pipe.", L"Pancake", MB_ICONERROR | MB_OK);
        return;
    }
    SetHandleInformation(read_pipe, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = write_pipe;
    startup.hStdError = write_pipe;
    startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);

    std::vector<wchar_t> mutable_command(command.begin(), command.end());
    mutable_command.push_back(L'\0');

    PROCESS_INFORMATION process_info{};
    const BOOL started = CreateProcessW(
        nullptr,
        mutable_command.data(),
        nullptr,
        nullptr,
        TRUE,
        CREATE_NO_WINDOW,
        nullptr,
        working_dir.c_str(),
        &startup,
        &process_info);

    CloseHandle(write_pipe);

    if (!started) {
        CloseHandle(read_pipe);
        MessageBoxW(g_hwnd, L"Could not start the solver process.", L"Pancake", MB_ICONERROR | MB_OK);
        return;
    }

    clear_log();
    append_log(L"> " + command + L"\r\n\r\n");
    close_solver_handles();
    g_solver.process_info = process_info;
    g_solver.stdout_read = read_pipe;
    g_solver.running = true;
    set_run_state(true);
    set_status(L"Solver running.");
    SetTimer(g_hwnd, kProcessTimer, 500, nullptr);
}

void stop_solver() {
    if (!g_solver.running) return;
    TerminateProcess(g_solver.process_info.hProcess, 1);
    append_log(L"\r\nSolver stop requested.\r\n");
    close_solver_handles();
    set_run_state(false);
    set_status(L"Solver stopped.");
}

bool parse_extent_values(const std::string& text, size_t extent_pos, int& x0, int& x1, int& y0, int& y1) {
    const size_t value_start = extent_pos + 8;
    const size_t value_end = text.find('"', value_start);
    if (value_end == std::string::npos) return false;

    std::stringstream values(text.substr(value_start, value_end - value_start));
    int z0 = 0;
    int z1 = 0;
    if (!(values >> x0 >> x1 >> y0 >> y1 >> z0 >> z1)) return false;
    return x1 > x0 && y1 > y0;
}

bool parse_vts_extents(const std::string& text, int& whole_nx, int& whole_ny, int& piece_x0, int& piece_x1, int& piece_y0, int& piece_y1) {
    int whole_x0 = 0;
    int whole_x1 = 0;
    int whole_y0 = 0;
    int whole_y1 = 0;
    const size_t whole_pos = text.find("WholeExtent=\"");
    if (whole_pos == std::string::npos || !parse_extent_values(text, whole_pos + 5, whole_x0, whole_x1, whole_y0, whole_y1)) {
        return false;
    }

    const size_t piece_tag = text.find("<Piece");
    if (piece_tag == std::string::npos) return false;
    const size_t piece_pos = text.find("Extent=\"", piece_tag);
    if (piece_pos == std::string::npos || !parse_extent_values(text, piece_pos, piece_x0, piece_x1, piece_y0, piece_y1)) {
        return false;
    }

    whole_nx = whole_x1 - whole_x0;
    whole_ny = whole_y1 - whole_y0;
    return whole_nx > 0 && whole_ny > 0;
}

bool parse_data_array(const std::string& text, const std::string& field, std::vector<double>& values) {
    const std::string name_token = "Name=\"" + field + "\"";
    const size_t name_pos = text.find(name_token);
    if (name_pos == std::string::npos) return false;

    const size_t tag_end = text.find('>', name_pos);
    if (tag_end == std::string::npos) return false;
    const size_t close = text.find("</DataArray>", tag_end);
    if (close == std::string::npos) return false;

    std::stringstream input(text.substr(tag_end + 1, close - tag_end - 1));
    values.clear();
    double value = 0.0;
    while (input >> value) values.push_back(value);
    return !values.empty();
}

void refresh_preview();

bool file_name_matches_vts(const fs::path& path, const std::string& prefix) {
    if (path.extension() != L".vts") return false;
    if (prefix.empty()) return true;
    const std::string name = wide_to_utf8(path.filename().wstring());
    return name.rfind(prefix + "_", 0) == 0;
}

bool parse_step_from_vts_name(const fs::path& path, const std::string& prefix, int& step) {
    const std::string name = wide_to_utf8(path.filename().wstring());
    const std::string front = prefix + "_";
    if (name.rfind(front, 0) != 0) return false;

    const size_t step_start = front.size();
    const size_t step_end = name.find('_', step_start);
    if (step_end == std::string::npos) return false;

    try {
        step = std::stoi(name.substr(step_start, step_end - step_start));
        return step >= 0;
    } catch (...) {
        return false;
    }
}

void collect_preview_steps_from_dir(const fs::path& directory, const std::string& prefix, std::vector<PreviewStep>& steps) {
    if (!fs::exists(directory)) return;

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        if (!file_name_matches_vts(entry.path(), prefix)) continue;
        int step = -1;
        if (!parse_step_from_vts_name(entry.path(), prefix, step)) continue;

        const auto existing = std::find_if(steps.begin(), steps.end(), [step](const PreviewStep& item) {
            return item.step == step;
        });
        if (existing == steps.end()) {
            steps.push_back({step, entry.path(), entry.last_write_time()});
        } else if (entry.last_write_time() > existing->write_time) {
            existing->path = entry.path();
            existing->write_time = entry.last_write_time();
        }
    }
}

void refresh_preview_step_items() {
    const int previous_step = [] {
        HWND combo = find_control(kComboPreviewStep);
        if (combo == nullptr) return -1;
        const LRESULT index = SendMessageW(combo, CB_GETCURSEL, 0, 0);
        if (index < 0) return -1;
        return static_cast<int>(SendMessageW(combo, CB_GETITEMDATA, static_cast<WPARAM>(index), 0));
    }();

    g_preview_steps.clear();
    const fs::path output_dir = output_directory_for_config(g_config);
    collect_preview_steps_from_dir(output_dir / L"flat" / L"processor0", g_config.filename_prefix, g_preview_steps);
    if (g_preview_steps.empty()) {
        collect_preview_steps_from_dir(output_dir / L"processor0", g_config.filename_prefix, g_preview_steps);
    }

    std::sort(g_preview_steps.begin(), g_preview_steps.end(), [](const PreviewStep& a, const PreviewStep& b) {
        return a.step < b.step;
    });

    HWND combo = find_control(kComboPreviewStep);
    if (combo == nullptr) return;

    g_loading = true;
    SendMessageW(combo, CB_RESETCONTENT, 0, 0);
    LRESULT latest_index = SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(L"Latest"));
    SendMessageW(combo, CB_SETITEMDATA, static_cast<WPARAM>(latest_index), static_cast<LPARAM>(-1));

    int selected_index = previous_step == -1 ? static_cast<int>(latest_index) : -1;
    for (const auto& item : g_preview_steps) {
        std::wstring label = L"step " + std::to_wstring(item.step);
        LRESULT index = SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(label.c_str()));
        SendMessageW(combo, CB_SETITEMDATA, static_cast<WPARAM>(index), static_cast<LPARAM>(item.step));
        if (item.step == previous_step) selected_index = static_cast<int>(index);
    }

    if (selected_index < 0) selected_index = static_cast<int>(latest_index);
    SendMessageW(combo, CB_SETCURSEL, static_cast<WPARAM>(selected_index), 0);
    g_loading = false;
}

int selected_preview_step() {
    HWND combo = find_control(kComboPreviewStep);
    if (combo == nullptr) return -1;
    const LRESULT index = SendMessageW(combo, CB_GETCURSEL, 0, 0);
    if (index < 0) return -1;
    return static_cast<int>(SendMessageW(combo, CB_GETITEMDATA, static_cast<WPARAM>(index), 0));
}

bool selected_preview_path(fs::path& path, int& step) {
    if (g_preview_steps.empty()) return false;

    step = selected_preview_step();
    if (step < 0) {
        const auto latest = std::max_element(g_preview_steps.begin(), g_preview_steps.end(), [](const PreviewStep& a, const PreviewStep& b) {
            return a.step < b.step;
        });
        if (latest == g_preview_steps.end()) return false;
        step = latest->step;
        path = latest->path;
        return true;
    }

    const auto item = std::find_if(g_preview_steps.begin(), g_preview_steps.end(), [step](const PreviewStep& candidate) {
        return candidate.step == step;
    });
    if (item == g_preview_steps.end()) return false;
    path = item->path;
    return true;
}

bool preview_file_for_step(const fs::path& path, const std::string& prefix, int step) {
    int file_step = -1;
    return file_name_matches_vts(path, prefix) && parse_step_from_vts_name(path, prefix, file_step) && file_step == step;
}

bool load_preview_values_for_step(
    const fs::path& seed_path,
    int step,
    const std::string& field,
    int& nx,
    int& ny,
    std::vector<double>& values,
    std::wstring& source_label) {
    const fs::path root = seed_path.parent_path().parent_path();
    if (!fs::exists(root)) return false;

    bool loaded_any = false;
    nx = 0;
    ny = 0;
    values.clear();
    std::vector<fs::path> step_files;

    for (const auto& entry : fs::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file()) continue;
        if (preview_file_for_step(entry.path(), g_config.filename_prefix, step)) {
            step_files.push_back(entry.path());
        }
    }
    std::sort(step_files.begin(), step_files.end());

    for (const auto& path : step_files) {
        const std::string text = read_file_text(path);
        int whole_nx = 0;
        int whole_ny = 0;
        int piece_x0 = 0;
        int piece_x1 = 0;
        int piece_y0 = 0;
        int piece_y1 = 0;
        std::vector<double> local_values;
        if (!parse_vts_extents(text, whole_nx, whole_ny, piece_x0, piece_x1, piece_y0, piece_y1) ||
            !parse_data_array(text, field, local_values)) {
            continue;
        }

        if (!loaded_any) {
            nx = whole_nx;
            ny = whole_ny;
            values.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), std::numeric_limits<double>::quiet_NaN());
            loaded_any = true;
        }
        if (whole_nx != nx || whole_ny != ny) continue;

        const int local_nx = piece_x1 - piece_x0;
        const int local_ny = piece_y1 - piece_y0;
        const size_t expected = static_cast<size_t>(local_nx) * static_cast<size_t>(local_ny);
        if (local_values.size() < expected) continue;

        for (int j = 0; j < local_ny; ++j) {
            const int global_j = piece_y0 + j;
            if (global_j < 0 || global_j >= ny) continue;
            for (int i = 0; i < local_nx; ++i) {
                const int global_i = piece_x0 + i;
                if (global_i < 0 || global_i >= nx) continue;
                values[static_cast<size_t>(global_j) * nx + global_i] =
                    local_values[static_cast<size_t>(j) * local_nx + i];
            }
        }
    }

    if (!loaded_any) return false;
    source_label = root.wstring() + L" (step " + std::to_wstring(step) + L")";
    return true;
}

void move_preview_step_selection(int delta) {
    HWND combo = find_control(kComboPreviewStep);
    if (combo == nullptr) return;
    const int count = static_cast<int>(SendMessageW(combo, CB_GETCOUNT, 0, 0));
    if (count <= 1) return;

    int index = static_cast<int>(SendMessageW(combo, CB_GETCURSEL, 0, 0));
    if (index <= 0) index = (delta < 0) ? count - 1 : 1;
    else index = std::clamp(index + delta, 1, count - 1);

    SendMessageW(combo, CB_SETCURSEL, static_cast<WPARAM>(index), 0);
    refresh_preview();
}

COLORREF heat_color(double t) {
    t = std::clamp(t, 0.0, 1.0);
    struct Stop { double at; int r; int g; int b; };
    constexpr std::array<Stop, 4> stops = {{
        {0.0, 33, 75, 99},
        {0.36, 35, 151, 139},
        {0.68, 241, 196, 83},
        {1.0, 188, 64, 56},
    }};

    Stop a = stops.front();
    Stop b = stops.back();
    for (size_t i = 1; i < stops.size(); ++i) {
        if (t <= stops[i].at) {
            a = stops[i - 1];
            b = stops[i];
            break;
        }
    }
    const double local = (b.at > a.at) ? (t - a.at) / (b.at - a.at) : 0.0;
    const int r = static_cast<int>(a.r + (b.r - a.r) * local);
    const int g = static_cast<int>(a.g + (b.g - a.g) * local);
    const int bl = static_cast<int>(a.b + (b.b - a.b) * local);
    return RGB(r, g, bl);
}

std::wstring compact_number_wide(double value) {
    return utf8_to_wide(compact_number(value));
}

std::string plot_number(double value) {
    if (!std::isfinite(value)) return "nan";
    const double magnitude = std::abs(value);
    std::ostringstream out;
    if (magnitude > 0.0 && (magnitude < 1.0e-3 || magnitude >= 1.0e5)) {
        out << std::setprecision(3) << std::scientific << value;
    } else {
        out << std::setprecision(5) << std::defaultfloat << value;
    }
    return out.str();
}

std::wstring plot_number_wide(double value) {
    return utf8_to_wide(plot_number(value));
}

std::wstring field_unit_label(const std::string& field) {
    if (field == "pressure") return L"Pa";
    if (field == "film_content") return L"-";
    if (field == "h") return L"m";
    if (field == "rho") return L"kg/m^3";
    if (field == "inlet_indicator") return L"-";
    if (field == "load_x" || field == "load_y" || field == "load_z") return L"N";
    if (field == "friction_torque") return L"N m";
    if (field == "velocity") return L"m/s";
    if (field == "theta_rad") return L"rad";
    if (field == "z_m") return L"m";
    return {};
}

std::wstring colorbar_title() {
    std::wstring title = utf8_to_wide(g_preview.field);
    const std::wstring unit = field_unit_label(g_preview.field);
    if (!unit.empty()) title += L" [" + unit + L"]";
    return title;
}

int preview_source_column(int display_column, int nx) {
    if (nx <= 0) return display_column;

    const double theta_zero_from_x_deg = 90.0 + g_config.load_angle_deg;
    const double cells = theta_zero_from_x_deg / 360.0 * static_cast<double>(nx);
    int shift = static_cast<int>(std::llround(cells));
    shift %= nx;
    if (shift < 0) shift += nx;

    return (display_column + shift) % nx;
}

void refresh_preview() {
    const std::wstring field_wide = get_combo_text(kComboPreviewField);
    const std::string field = normalise_key(wide_to_utf8(field_wide));
    if (field.empty()) return;

    refresh_preview_step_items();

    fs::path selected_path;
    int selected_step = -1;
    if (!selected_preview_path(selected_path, selected_step)) {
        g_preview = PreviewData{};
        g_preview.message = L"No VTS results found yet.";
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }

    int nx = 0;
    int ny = 0;
    std::vector<double> values;
    std::wstring source_label;
    if (!load_preview_values_for_step(selected_path, selected_step, field, nx, ny, values, source_label)) {
        g_preview = PreviewData{};
        g_preview.message = L"Selected field was not found in the selected timestep.";
        g_preview.source_path = selected_path.wstring();
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }

    const size_t expected = static_cast<size_t>(nx) * static_cast<size_t>(ny);
    if (values.size() < expected) {
        g_preview = PreviewData{};
        g_preview.message = L"Latest VTS file has incomplete field data.";
        g_preview.source_path = selected_path.wstring();
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }
    values.resize(expected);

    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (!std::isfinite(value)) continue;
        min_value = std::min(min_value, value);
        max_value = std::max(max_value, value);
    }
    if (!std::isfinite(min_value) || !std::isfinite(max_value)) {
        min_value = 0.0;
        max_value = 0.0;
    }

    g_preview.ok = true;
    g_preview.nx = nx;
    g_preview.ny = ny;
    g_preview.step = selected_step;
    g_preview.min_value = min_value;
    g_preview.max_value = max_value;
    g_preview.field = field;
    g_preview.values = std::move(values);
    g_preview.source_path = source_label;
    InvalidateRect(g_preview_canvas, nullptr, TRUE);
}

void poll_solver() {
    if (!g_solver.running) return;

    char buffer[4096];
    DWORD available = 0;
    while (g_solver.stdout_read != nullptr &&
           PeekNamedPipe(g_solver.stdout_read, nullptr, 0, nullptr, &available, nullptr) &&
           available > 0) {
        DWORD bytes_read = 0;
        const DWORD to_read = std::min<DWORD>(available, sizeof(buffer) - 1);
        if (!ReadFile(g_solver.stdout_read, buffer, to_read, &bytes_read, nullptr) || bytes_read == 0) break;
        buffer[bytes_read] = '\0';
        append_log(utf8_to_wide(std::string(buffer, bytes_read)));
    }

    refresh_preview();

    DWORD exit_code = STILL_ACTIVE;
    if (GetExitCodeProcess(g_solver.process_info.hProcess, &exit_code) && exit_code != STILL_ACTIVE) {
        append_log(L"\r\nSolver exited with code " + std::to_wstring(exit_code) + L".\r\n");
        close_solver_handles();
        set_run_state(false);
        set_status(exit_code == 0 ? L"Solver finished." : L"Solver exited with errors.");
        KillTimer(g_hwnd, kProcessTimer);
        refresh_preview();
    }
}

void open_results_folder() {
    if (!save_config()) return;
    const fs::path output_dir = output_directory_for_config(g_config);
    fs::create_directories(output_dir);
    ShellExecuteW(g_hwnd, L"open", output_dir.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
}

void open_paraview_result() {
    if (!save_config()) return;
    const fs::path output_dir = output_directory_for_config(g_config);
    fs::path pvd;
    if (g_config.output_write_flat && fs::exists(output_dir / L"flat" / L"results.pvd")) {
        pvd = output_dir / L"flat" / L"results.pvd";
    } else {
        pvd = output_dir / L"results.pvd";
    }

    if (!fs::exists(pvd)) {
        MessageBoxW(g_hwnd, L"No PVD result file exists yet.", L"Pancake", MB_ICONINFORMATION | MB_OK);
        return;
    }
    ShellExecuteW(g_hwnd, L"open", pvd.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
}

void sync_raw_from_form() {
    SimulationConfig cfg = g_config;
    if (!apply_form_to_config(cfg)) return;
    g_config = cfg;
    g_loading = true;
    set_edit_text(kEditRawConfig, form_config_text(g_config));
    g_loading = false;
    set_status(L"Raw config synchronized from the form.");
}

void apply_raw_to_form() {
    const std::string raw_text = get_edit_utf8(kEditRawConfig);
    SimulationConfig cfg;
    try {
        cfg.load_from_text(raw_text);
    } catch (const std::exception& ex) {
        MessageBoxW(g_hwnd, utf8_to_wide(ex.what()).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return;
    }
    g_config = cfg;
    g_config_value_text = parse_config_value_text(raw_text);
    g_config_inlets_text = parse_config_inlets_text(raw_text);
    populate_form_from_config(g_config);
    g_dirty = true;
    set_status(L"Raw config applied to the form.");
}

void layout_controls();

void show_tab(int tab) {
    g_active_tab = tab;
    for (int i = 0; i < kTabCount; ++i) {
        for (HWND child : g_tab_controls[i]) {
            ShowWindow(child, i == tab ? SW_SHOW : SW_HIDE);
        }
    }
    if (tab == kRawTab) {
        sync_raw_from_form();
    }
    if (tab == kRunTab) {
        refresh_preview();
    }
    raise_tab_page_controls(tab);
}

RECT run_content_rect(int width, int height) {
    RECT rect{};
    rect.left = 32;
    rect.top = 118;
    rect.right = static_cast<LONG>(std::max(static_cast<int>(rect.left) + 720, width - 32));
    rect.bottom = static_cast<LONG>(std::max(static_cast<int>(rect.top) + 430, height - 28));
    return rect;
}

int clamp_run_bottom_height(int requested_height, int content_height) {
    constexpr int min_bottom = 150;
    constexpr int min_view = 260;
    const int max_bottom = std::max(min_bottom, content_height - min_view - kSplitterWidth - kPaneGap);
    return std::clamp(requested_height, min_bottom, max_bottom);
}

bool param_item_visible(const ParamItem& item) {
    if (item.visibility == ParamVisibility::CircularInletOnly) {
        return combo_is(kComboInletType, L"CIRCULAR");
    }
    return true;
}

void show_param_item(const ParamItem& item, bool visible) {
    const int command = visible ? SW_SHOW : SW_HIDE;
    if (item.group != nullptr) ShowWindow(item.group, command);
    if (item.label != nullptr) ShowWindow(item.label, command);
    if (item.value != nullptr) ShowWindow(item.value, command);
    if (item.unit != nullptr) ShowWindow(item.unit, command);
}

void layout_parameter_panel() {
    if (g_param_panel == nullptr) return;

    RECT rect{};
    GetClientRect(g_param_panel, &rect);
    const int panel_width = std::max(260, static_cast<int>(rect.right - rect.left));
    g_param_panel_width = panel_width;

    const int pad = 12;
    const int group_x = 8;
    const int group_width = std::max(120, panel_width - 24);
    const int label_x = 22;
    const int unit_width = 74;
    const int gap = 8;
    const int value_width = std::max(72, group_width - 28 - 116 - unit_width - gap);
    const int value_x = label_x + 112;
    const int unit_x = value_x + value_width + gap;
    const int row_height = 30;

    int y = pad - g_param_scroll_y;
    ParamItem* open_group = nullptr;
    int open_group_top = 0;

    auto finish_group = [&](int bottom_y) {
        if (open_group == nullptr) return;
        const int height = std::max(44, bottom_y - open_group_top + 8);
        move_control(open_group->group, group_x, open_group_top, group_width, height);
    };

    for (ParamItem& item : g_param_items) {
        if (item.kind == ParamItemKind::Group) {
            finish_group(y);
            open_group = &item;
            open_group_top = y;
            show_param_item(item, true);
            y += 30;
            continue;
        }

        const bool visible = param_item_visible(item);
        show_param_item(item, visible);
        if (!visible) continue;

        if (item.kind == ParamItemKind::CheckRow) {
            move_control(item.value, label_x, y, group_width - 36, 24);
            y += row_height;
            continue;
        }

        move_control(item.label, label_x, y + 4, 110, 22);
        move_control(item.value, value_x, y, value_width, 24);
        if (item.kind == ParamItemKind::ComboRow) {
            move_control(item.value, value_x, y - 1, value_width + unit_width + gap, 220);
        } else {
            move_control(item.unit, unit_x, y - 1, unit_width, 220);
        }
        y += row_height;
    }
    finish_group(y);

    g_param_content_height = y + g_param_scroll_y + pad;
    const int visible_height = rect.bottom - rect.top;
    const int max_scroll = std::max(0, g_param_content_height - visible_height);
    g_param_scroll_y = std::clamp(g_param_scroll_y, 0, max_scroll);

    SCROLLINFO info{};
    info.cbSize = sizeof(info);
    info.fMask = SIF_RANGE | SIF_PAGE | SIF_POS;
    info.nMin = 0;
    info.nMax = std::max(g_param_content_height, visible_height);
    info.nPage = static_cast<UINT>(std::max(1, visible_height));
    info.nPos = g_param_scroll_y;
    SetScrollInfo(g_param_panel, SB_VERT, &info, TRUE);
}

void scroll_parameter_panel(int delta) {
    if (g_param_panel == nullptr) return;
    RECT rect{};
    GetClientRect(g_param_panel, &rect);
    const int max_scroll = std::max(0, g_param_content_height - static_cast<int>(rect.bottom - rect.top));
    const int next = std::clamp(g_param_scroll_y + delta, 0, max_scroll);
    if (next == g_param_scroll_y) return;
    g_param_scroll_y = next;
    layout_parameter_panel();
    InvalidateRect(g_param_panel, nullptr, TRUE);
}

void update_run_splitter_from_client_y(int client_y) {
    RECT client{};
    GetClientRect(g_hwnd, &client);
    const RECT content = run_content_rect(client.right - client.left, client.bottom - client.top);
    const int content_height = content.bottom - content.top;
    const int requested_height = content.bottom - client_y - kSplitterWidth - kPaneGap;
    g_run_bottom_height = clamp_run_bottom_height(requested_height, content_height);
    layout_controls();
}

void layout_run_tab(const RECT& content) {
    const int content_width = content.right - content.left;
    const int content_height = content.bottom - content.top;
    if (content_width <= 0 || content_height <= 0) return;

    const int parameter_width = std::clamp(static_cast<int>(content_width * 0.32), 340, 440);
    const int left_width = std::max(420, content_width - parameter_width - kPaneGap);
    const int parameter_left = content.left + left_width + kPaneGap;
    move_control(g_param_panel, parameter_left, content.top, parameter_width, content_height);
    layout_parameter_panel();

    g_run_bottom_height = clamp_run_bottom_height(g_run_bottom_height, content_height);

    const int result_height = content_height - g_run_bottom_height - kSplitterWidth - kPaneGap;
    const int splitter_top = content.top + result_height + kPaneGap;
    const int bottom_top = splitter_top + kSplitterWidth + kPaneGap;

    move_control(g_run_result_group, content.left, content.top, left_width, result_height);
    move_control(g_run_field_label, content.left + 18, content.top + 32, 42, 22);
    move_control_id(kComboPreviewField, content.left + 64, content.top + 28, 160, 220);
    move_control(g_run_step_label, content.left + 238, content.top + 32, 42, 22);
    move_control_id(kComboPreviewStep, content.left + 284, content.top + 28, 110, 220);
    move_control_id(kPrevStepId, content.left + 404, content.top + 26, 34, 30);
    move_control_id(kNextStepId, content.left + 444, content.top + 26, 34, 30);
    move_control_id(kRefreshPreviewId, content.left + 490, content.top + 26, 86, 30);
    move_control_id(kOpenParaviewId, content.left + 586, content.top + 26, 90, 30);

    const int canvas_top = content.top + 68;
    const int canvas_height = content.top + result_height - canvas_top - 20;
    move_control(g_preview_canvas, content.left + 18, canvas_top, left_width - 36, canvas_height);

    move_control(g_run_splitter, content.left, splitter_top, left_width, kSplitterWidth);

    move_control(g_run_log_group, content.left, bottom_top, left_width, g_run_bottom_height);
    move_control_id(kEditLog, content.left + 18, bottom_top + 28, left_width - 36, g_run_bottom_height - 46);
    bring_control_id_to_front(kComboPreviewField);
    bring_control_id_to_front(kComboPreviewStep);
    bring_control_id_to_front(kPrevStepId);
    bring_control_id_to_front(kNextStepId);
    bring_control_id_to_front(kRefreshPreviewId);
    bring_control_id_to_front(kOpenParaviewId);
}

void layout_controls() {
    RECT client{};
    GetClientRect(g_hwnd, &client);
    const int width = client.right - client.left;
    const int height = client.bottom - client.top;

    const int margin = 14;
    const int toolbar_y = 14;
    int x = margin;
    move_control_id(kSaveId, x, toolbar_y, 92, 30);
    x += 98;
    move_control_id(kReloadId, x, toolbar_y, 92, 30);
    x += 98;
    move_control_id(kDefaultsId, x, toolbar_y, 92, 30);
    x += 106;
    move_control_id(kOpenResultsId, x, toolbar_y, 118, 30);

    const int solver_width = std::min(470, std::max(410, width - 500));
    const int solver_left = std::max(500, width - solver_width - margin);
    move_control(g_run_solver_group, solver_left, 8, solver_width, 64);
    move_control(g_run_solver_label, solver_left + 14, 34, 28, 22);
    move_control_id(kEditSolverPath, solver_left + 46, 30, solver_width - 286, 24);
    move_control(g_run_ranks_label, solver_left + solver_width - 234, 34, 46, 22);
    move_control_id(kEditMpiRanks, solver_left + solver_width - 184, 30, 44, 24);
    move_control(g_run_button, solver_left + solver_width - 130, 28, 58, 28);
    move_control(g_stop_button, solver_left + solver_width - 66, 28, 54, 28);
    bring_control_to_front(g_run_button);
    bring_control_to_front(g_stop_button);

    move_control(g_status, margin, 48, std::max(240, solver_left - margin - 12), 22);

    move_control(g_tabs, margin, 78, width - 2 * margin, height - 92);

    const int content_top = 118;
    const int content_left = 32;
    const int content_right = width - 32;
    const int content_bottom = height - 28;

    HWND raw = find_control(kEditRawConfig);
    move_control(raw, content_left, content_top + 42, content_right - content_left, content_bottom - content_top - 52);
    move_control_id(kSyncRawId, content_left, content_top + 4, 128, 30);
    move_control_id(kApplyRawId, content_left + 138, content_top + 4, 128, 30);

    layout_run_tab(run_content_rect(width, height));
    raise_tab_page_controls(g_active_tab);
}

void create_parameter_tab() {
    g_param_panel = create_window_control(
        kWorkspaceTab,
        kParamPanelClass,
        L"",
        WS_VSCROLL | WS_CLIPCHILDREN,
        WS_EX_CLIENTEDGE,
        0,
        0,
        1,
        1,
        0);

    g_parent_override = g_param_panel;
    g_param_items.clear();

    param_group(L"Geometry");
    param_unit_row(L"R", kEditR, kComboUnitR, {L"m", L"cm", L"mm", L"um"});
    param_unit_row(L"c", kEditC, kComboUnitC, {L"m", L"cm", L"mm", L"um"});
    param_unit_row(L"e", kEditE, kComboUnitE, {L"m", L"cm", L"mm", L"um"});
    param_unit_row(L"L", kEditL, kComboUnitL, {L"m", L"cm", L"mm", L"um"});
    param_unit_row(L"attitude", kEditAttitude, kComboUnitAttitude, {L"deg", L"rad"});
    param_unit_row(L"load angle", kEditLoadAngle, kComboUnitLoadAngle, {L"deg", L"rad"});
    param_unit_row(L"tilt_x", kEditTiltX, kComboUnitTiltX, {L"m/m", L"deg", L"rad"});
    param_unit_row(L"tilt_y", kEditTiltY, kComboUnitTiltY, {L"m/m", L"deg", L"rad"});

    param_group(L"Grid and Time");
    param_unit_row(L"n_theta", kEditNTheta, kComboUnitNTheta, {L"cells"});
    param_unit_row(L"n_z", kEditNZ, kComboUnitNZ, {L"cells"});
    param_unit_row(L"end_t", kEditEndT, kComboUnitEndT, {L"s", L"ms"});
    param_unit_row(L"dt", kEditDt, kComboUnitDt, {L"s", L"ms"});
    param_unit_row(L"write interval", kEditWriteInterval, kComboUnitWriteInterval, {L"s", L"ms"});

    param_group(L"Physics");
    param_unit_row(L"omega", kEditOmega, kComboUnitOmega, {L"rad/s", L"rpm"});
    param_unit_row(L"mu", kEditMu, kComboUnitMu, {L"Pa s"});
    param_unit_row(L"rho", kEditRho, kComboUnitRho, {L"kg/m3"});
    param_unit_row(L"p_cav", kEditPCav, kComboUnitPCav, {L"Pa", L"kPa", L"MPa"});
    param_unit_row(L"bulk modulus", kEditBulkModulus, kComboUnitBulkModulus, {L"Pa", L"kPa", L"MPa"});

    param_group(L"Cavitation");
    param_combo_row(L"model", kComboCavitation, {L"ELROD_ADAMS", L"GUMBEL", L"FULL_SOMMERFELD"});
    param_unit_row(L"max outer", kEditMaxOuter, kComboUnitMaxOuter, {L"iter"});
    param_unit_row(L"outer tol", kEditOuterTol, kComboUnitOuterTol, {L"-"});
    param_unit_row(L"theta min", kEditThetaMin, kComboUnitThetaMin, {L"-"});
    param_check_row(kCheckLogOuter, L"log outer iterations");

    param_group(L"South Boundary");
    param_combo_row(L"type", kComboBcSouth, {L"DIRICHLET", L"NEUMANN", L"INLET_OUTLET"});
    param_unit_row(L"pressure", kEditBcSouthVal, kComboUnitBcSouthVal, {L"Pa", L"kPa", L"MPa"});
    param_unit_row(L"film content", kEditBcSouthTheta, kComboUnitBcSouthTheta, {L"-"});

    param_group(L"North Boundary");
    param_combo_row(L"type", kComboBcNorth, {L"DIRICHLET", L"NEUMANN", L"INLET_OUTLET"});
    param_unit_row(L"pressure", kEditBcNorthVal, kComboUnitBcNorthVal, {L"Pa", L"kPa", L"MPa"});
    param_unit_row(L"film content", kEditBcNorthTheta, kComboUnitBcNorthTheta, {L"-"});

    param_group(L"Inlet");
    param_combo_row(L"type", kComboInletType, {L"GROOVE", L"CIRCULAR"});
    param_unit_row(L"theta", kEditInletTheta, kComboUnitInletTheta, {L"deg", L"rad"});
    param_unit_row(L"z", kEditInletZ, kComboUnitInletZ, {L"m", L"cm", L"mm", L"um"}, ParamVisibility::CircularInletOnly);
    g_inlet_z_label = g_param_items.back().label;
    param_unit_row(L"width", kEditInletSize, kComboUnitInletSize, {L"deg", L"rad"});
    g_inlet_size_label = g_param_items.back().label;
    param_unit_row(L"p_supply", kEditInletPressure, kComboUnitInletPressure, {L"Pa", L"kPa", L"MPa"});

    g_parent_override = nullptr;
}

void create_output_tab() {
    const int top = 118;
    create_group(kOutputTab, L"Files", 32, top, 520, 150);
    create_label(kOutputTab, L"output_dir", 52, top + 34, 126);
    create_edit(kOutputTab, kEditOutputDir, 178, top + 30, 320);
    create_label(kOutputTab, L"filename_prefix", 52, top + 68, 126);
    create_edit(kOutputTab, kEditFilenamePrefix, 178, top + 64, 320);
    create_checkbox(kOutputTab, kCheckOutput3d, L"write curved 3D VTK/PVD", 178, top + 100, 210);
    create_checkbox(kOutputTab, kCheckOutputFlat, L"write flat unwrapped VTK/PVD", 178, top + 124, 230);

    create_group(kOutputTab, L"Data to Keep", 582, top, 440, 252);
    int field_x = 604;
    int field_y = top + 32;
    for (size_t i = 0; i < kOutputFields.size(); ++i) {
        const int col = static_cast<int>(i / 5);
        const int row = static_cast<int>(i % 5);
        create_checkbox(kOutputTab, kOutputFields[i].id, kOutputFields[i].label, field_x + col * 210, field_y + row * 32, 190);
    }
}

void create_run_tab() {
    const int top = 118;
    g_run_result_group = create_group(kRunTab, L"Results", 32, top, 720, 372);
    g_run_field_label = create_label(kRunTab, L"field", 52, top + 32, 42);
    create_combo(
        kRunTab,
        kComboPreviewField,
        98,
        top + 28,
        180,
        {L"pressure", L"film_content", L"h", L"rho", L"inlet_indicator", L"load_x", L"load_y", L"load_z", L"friction_torque"});
    g_run_step_label = create_label(kRunTab, L"step", 292, top + 32, 42);
    create_combo(kRunTab, kComboPreviewStep, 338, top + 28, 126, {L"Latest"});
    create_tab_button(kRunTab, kPrevStepId, L"<", 474, top + 26, 34);
    create_tab_button(kRunTab, kNextStepId, L">", 514, top + 26, 34);
    create_tab_button(kRunTab, kRefreshPreviewId, L"Refresh", 560, top + 26, 92);
    create_tab_button(kRunTab, kOpenParaviewId, L"Open PVD", 662, top + 26, 96);

    g_preview_canvas = create_window_control(
        kRunTab,
        kPreviewClass,
        L"",
        0,
        WS_EX_CLIENTEDGE,
        52,
        top + 68,
        680,
        284,
        kPreviewCanvasId);

    g_run_splitter = create_window_control(kRunTab, kSplitterClass, L"", 0, 0, 756, top, kSplitterWidth, 372, 0);

    g_run_log_group = create_group(kRunTab, L"Process Log", 772, top + 190, 250, 352);
    create_window_control(
        kRunTab,
        g_richedit_available ? L"RICHEDIT50W" : L"EDIT",
        L"",
        WS_TABSTOP | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY | WS_VSCROLL,
        WS_EX_CLIENTEDGE,
        52,
        top + 420,
        950,
        104,
        kEditLog);
    HWND log = find_control(kEditLog);
    if (log != nullptr) {
        SendMessageW(log, WM_SETFONT, reinterpret_cast<WPARAM>(g_mono_font), TRUE);
        if (g_richedit_available) {
            SendMessageW(log, EM_SETBKGNDCOLOR, 0, RGB(250, 251, 253));
        }
    }
}

void create_raw_tab() {
    const int top = 118;
    create_tab_button(kRawTab, kSyncRawId, L"Sync From Form", 32, top + 4, 128);
    create_tab_button(kRawTab, kApplyRawId, L"Apply Raw", 170, top + 4, 128);
    create_edit(kRawTab, kEditRawConfig, 32, top + 46, 960, 430, ES_MULTILINE | ES_AUTOVSCROLL | ES_AUTOHSCROLL | ES_WANTRETURN | WS_VSCROLL | WS_HSCROLL);
    HWND raw = find_control(kEditRawConfig);
    if (raw != nullptr) SendMessageW(raw, WM_SETFONT, reinterpret_cast<WPARAM>(g_mono_font), TRUE);
}

LRESULT CALLBACK tab_subclass_proc(HWND, UINT, WPARAM, LPARAM, UINT_PTR, DWORD_PTR);

void create_controls() {
    create_button(kSaveId, L"Save", 14, 14, 92);
    create_button(kReloadId, L"Reload", 112, 14, 92);
    create_button(kDefaultsId, L"Defaults", 210, 14, 92);
    create_button(kOpenResultsId, L"Open Results", 316, 14, 118);

    g_status = CreateWindowExW(
        0,
        L"STATIC",
        L"",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | SS_LEFT,
        14,
        48,
        860,
        22,
        g_hwnd,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kStatusId)),
        GetModuleHandleW(nullptr),
        nullptr);
    SendMessageW(g_status, WM_SETFONT, reinterpret_cast<WPARAM>(g_font), TRUE);
    remember_control(kStatusId, g_status);

    g_run_solver_group = create_group(-1, L"Solver Backend", 650, 8, 430, 64);
    g_run_solver_label = create_label(-1, L"exe", 664, 34, 34);
    create_edit(-1, kEditSolverPath, 700, 30, 184);
    g_run_ranks_label = create_label(-1, L"ranks", 894, 34, 44);
    create_edit(-1, kEditMpiRanks, 940, 30, 42);
    g_run_button = create_button(kRunId, L"Run", 992, 28, 40, 28);
    g_stop_button = create_button(kStopId, L"Stop", 1038, 28, 42, 28);
    EnableWindow(g_stop_button, FALSE);

    g_tabs = CreateWindowExW(
        0,
        WC_TABCONTROLW,
        L"",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_TABSTOP,
        14,
        78,
        1060,
        590,
        g_hwnd,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kTabsId)),
        GetModuleHandleW(nullptr),
        nullptr);
    SendMessageW(g_tabs, WM_SETFONT, reinterpret_cast<WPARAM>(g_font), TRUE);
    SetWindowSubclass(g_tabs, tab_subclass_proc, 0, 0);

    TCITEMW item{};
    item.mask = TCIF_TEXT;
    item.pszText = const_cast<wchar_t*>(L"Workspace");
    TabCtrl_InsertItem(g_tabs, kWorkspaceTab, &item);
    item.pszText = const_cast<wchar_t*>(L"Output");
    TabCtrl_InsertItem(g_tabs, kOutputTab, &item);
    item.pszText = const_cast<wchar_t*>(L"Raw Config");
    TabCtrl_InsertItem(g_tabs, kRawTab, &item);

    create_parameter_tab();
    create_output_tab();
    create_run_tab();
    create_raw_tab();
    populate_preview_fields();

    // Group-box frames are created interleaved with the combos they enclose;
    // force every combo above all frames so the closed control always paints
    // (previously some only appeared after a click).
    for (const auto& ref : g_control_refs) {
        wchar_t class_name[16] = {};
        GetClassNameW(ref.hwnd, class_name, 15);
        if (lstrcmpiW(class_name, L"COMBOBOX") == 0) bring_control_to_front(ref.hwnd);
    }

    TabCtrl_SetCurSel(g_tabs, kWorkspaceTab);
    show_tab(kWorkspaceTab);
}

void create_fonts() {
    g_font = CreateFontW(
        -16, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
        DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    g_heading_font = CreateFontW(
        -16, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
        DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    g_mono_font = CreateFontW(
        -15, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
        DEFAULT_PITCH | FF_MODERN, L"Consolas");
}

void draw_preview(HWND hwnd, HDC hdc) {
    RECT rect{};
    GetClientRect(hwnd, &rect);
    FillRect(hdc, &rect, g_brush_field);

    SetBkMode(hdc, TRANSPARENT);
    SelectObject(hdc, g_font);
    SetTextColor(hdc, kColText);

    if (!g_preview.ok) {
        DrawTextW(hdc, g_preview.message.c_str(), -1, &rect, DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
        return;
    }

    // Reserve generous margins so the title/stat row and axis labels never
    // sit on top of the contour.
    RECT plot = rect;
    plot.left += 64;
    plot.top += 82;
    plot.right -= 178;
    plot.bottom -= 78;
    if (plot.right <= plot.left || plot.bottom <= plot.top) return;

    std::wstring title = utf8_to_wide(g_preview.field) + L"  \x2022  step " + std::to_wstring(g_preview.step);
    RECT title_rect{rect.left + 16, rect.top + 12, rect.right - 16, rect.top + 34};
    SelectObject(hdc, g_heading_font);
    SetTextColor(hdc, kColText);
    DrawTextW(hdc, title.c_str(), -1, &title_rect, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    SelectObject(hdc, g_font);

    std::ostringstream stats;
    stats << "min " << plot_number(g_preview.min_value)
          << "      max " << plot_number(g_preview.max_value)
          << "      grid " << g_preview.nx << " x " << g_preview.ny;
    RECT stats_rect{rect.left + 16, rect.top + 42, rect.right - 16, rect.top + 62};
    SetTextColor(hdc, kColMuted);
    DrawTextW(hdc, utf8_to_wide(stats.str()).c_str(), -1, &stats_rect, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);

    const double range = g_preview.max_value - g_preview.min_value;
    for (int j = 0; j < g_preview.ny; ++j) {
        for (int i = 0; i < g_preview.nx; ++i) {
            const int source_i = preview_source_column(i, g_preview.nx);
            const double value = g_preview.values[static_cast<size_t>(j) * g_preview.nx + source_i];
            const double t = (range > 0.0 && std::isfinite(value)) ? (value - g_preview.min_value) / range : 0.5;
            HBRUSH brush = CreateSolidBrush(heat_color(t));
            RECT cell{};
            cell.left = plot.left + (i * (plot.right - plot.left)) / g_preview.nx;
            cell.right = plot.left + ((i + 1) * (plot.right - plot.left)) / g_preview.nx + 1;
            const int flipped_j = g_preview.ny - 1 - j;
            cell.top = plot.top + (flipped_j * (plot.bottom - plot.top)) / g_preview.ny;
            cell.bottom = plot.top + ((flipped_j + 1) * (plot.bottom - plot.top)) / g_preview.ny + 1;
            FillRect(hdc, &cell, brush);
            DeleteObject(brush);
        }
    }

    HPEN border = CreatePen(PS_SOLID, 1, kColBorder);
    HGDIOBJ old_pen = SelectObject(hdc, border);
    HGDIOBJ old_brush = SelectObject(hdc, GetStockObject(NULL_BRUSH));
    Rectangle(hdc, plot.left, plot.top, plot.right, plot.bottom);
    SelectObject(hdc, old_brush);
    SelectObject(hdc, old_pen);
    DeleteObject(border);

    HPEN axis_pen = CreatePen(PS_SOLID, 1, RGB(120, 130, 143));
    old_pen = SelectObject(hdc, axis_pen);
    MoveToEx(hdc, plot.left, plot.bottom, nullptr);
    LineTo(hdc, plot.right, plot.bottom);
    MoveToEx(hdc, plot.left, plot.top, nullptr);
    LineTo(hdc, plot.left, plot.bottom);
    SetTextColor(hdc, kColMuted);
    for (int tick = 0; tick <= 4; ++tick) {
        const int x = plot.left + tick * (plot.right - plot.left) / 4;
        MoveToEx(hdc, x, plot.bottom, nullptr);
        LineTo(hdc, x, plot.bottom + 5);
        const std::wstring label = std::to_wstring(tick * 90);
        RECT tick_rect{x - 24, plot.bottom + 6, x + 24, plot.bottom + 24};
        DrawTextW(hdc, label.c_str(), -1, &tick_rect, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
    }
    for (int tick = 0; tick <= 4; ++tick) {
        const int y = plot.bottom - tick * (plot.bottom - plot.top) / 4;
        MoveToEx(hdc, plot.left - 5, y, nullptr);
        LineTo(hdc, plot.left, y);
        const double z_value = g_config.L * static_cast<double>(tick) / 4.0;
        const std::wstring label = plot_number_wide(z_value);
        RECT tick_rect{rect.left + 4, y - 10, plot.left - 8, y + 10};
        DrawTextW(hdc, label.c_str(), -1, &tick_rect, DT_RIGHT | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
    }
    RECT x_label{plot.left, plot.bottom + 26, plot.right, plot.bottom + 44};
    DrawTextW(hdc, L"angle from +Y load reference (deg)", -1, &x_label, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
    HFONT y_font = CreateFontW(
        -16, 0, 900, 900, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
        DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    HGDIOBJ old_y_font = SelectObject(hdc, y_font);
    const UINT old_align = SetTextAlign(hdc, TA_CENTER | TA_BASELINE);
    TextOutW(hdc, rect.left + 18, (plot.top + plot.bottom) / 2, L"z [m]", 5);
    SetTextAlign(hdc, old_align);
    SelectObject(hdc, old_y_font);
    DeleteObject(y_font);
    SelectObject(hdc, old_pen);
    DeleteObject(axis_pen);

    RECT colorbar{};
    colorbar.left = plot.right + 24;
    colorbar.right = colorbar.left + 18;
    colorbar.top = plot.top;
    colorbar.bottom = plot.bottom;
    const int colorbar_height = std::max(1, static_cast<int>(colorbar.bottom - colorbar.top));
    for (int y = colorbar.top; y < colorbar.bottom; ++y) {
        const double t = 1.0 - static_cast<double>(y - colorbar.top) / colorbar_height;
        HPEN color_pen = CreatePen(PS_SOLID, 1, heat_color(t));
        HGDIOBJ old_color_pen = SelectObject(hdc, color_pen);
        MoveToEx(hdc, colorbar.left, y, nullptr);
        LineTo(hdc, colorbar.right, y);
        SelectObject(hdc, old_color_pen);
        DeleteObject(color_pen);
    }
    HPEN colorbar_border = CreatePen(PS_SOLID, 1, kColBorder);
    HGDIOBJ old_cb_pen = SelectObject(hdc, colorbar_border);
    HGDIOBJ old_color_brush = SelectObject(hdc, GetStockObject(NULL_BRUSH));
    Rectangle(hdc, colorbar.left, colorbar.top, colorbar.right, colorbar.bottom);
    SelectObject(hdc, old_color_brush);
    SelectObject(hdc, old_cb_pen);
    DeleteObject(colorbar_border);

    // Colorbar title and tick labels sit to the right with enough room for
    // scientific notation.
    RECT color_title{plot.right + 8, colorbar.top - 26, rect.right - 8, colorbar.top - 6};
    SelectObject(hdc, g_heading_font);
    SetTextColor(hdc, kColText);
    const std::wstring cb_title = colorbar_title();
    DrawTextW(hdc, cb_title.c_str(), -1, &color_title, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    SelectObject(hdc, g_font);
    HPEN tick_pen = CreatePen(PS_SOLID, 1, RGB(120, 130, 143));
    HGDIOBJ old_tick_pen = SelectObject(hdc, tick_pen);
    HFONT color_tick_font = CreateFontW(
        -13, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
        DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    HGDIOBJ old_tick_font = SelectObject(hdc, color_tick_font);
    const int n_color_ticks = 9;
    for (int tick = 0; tick < n_color_ticks; ++tick) {
        const double frac = (n_color_ticks == 1) ? 0.0 : static_cast<double>(tick) / (n_color_ticks - 1);
        const int y = colorbar.bottom - static_cast<int>(std::llround(frac * (colorbar.bottom - colorbar.top)));
        const double value = g_preview.min_value + frac * (g_preview.max_value - g_preview.min_value);
        MoveToEx(hdc, colorbar.right, y, nullptr);
        LineTo(hdc, colorbar.right + 5, y);
        RECT tick_label{colorbar.right + 8, y - 9, rect.right - 8, y + 10};
        DrawTextW(hdc, plot_number_wide(value).c_str(), -1, &tick_label, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
    }
    SelectObject(hdc, old_tick_font);
    DeleteObject(color_tick_font);
    SelectObject(hdc, old_tick_pen);
    DeleteObject(tick_pen);

    RECT source_rect{rect.left + 16, plot.bottom + 50, rect.right - 16, plot.bottom + 70};
    SetTextColor(hdc, kColMuted);
    DrawTextW(hdc, g_preview.source_path.c_str(), -1, &source_rect, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
}

LRESULT CALLBACK preview_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_PAINT: {
            PAINTSTRUCT paint{};
            HDC hdc = BeginPaint(hwnd, &paint);
            draw_preview(hwnd, hdc);
            EndPaint(hwnd, &paint);
            return 0;
        }
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK splitter_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_SETCURSOR:
            SetCursor(LoadCursorW(nullptr, IDC_SIZENS));
            return TRUE;

        case WM_LBUTTONDOWN:
            g_run_splitter_dragging = true;
            SetCapture(hwnd);
            return 0;

        case WM_MOUSEMOVE:
            if (g_run_splitter_dragging) {
                POINT point{};
                GetCursorPos(&point);
                ScreenToClient(g_hwnd, &point);
                update_run_splitter_from_client_y(point.y);
                return 0;
            }
            break;

        case WM_LBUTTONUP:
            if (g_run_splitter_dragging) {
                g_run_splitter_dragging = false;
                ReleaseCapture();
                return 0;
            }
            break;

        case WM_CAPTURECHANGED:
            g_run_splitter_dragging = false;
            return 0;

        case WM_PAINT: {
            PAINTSTRUCT paint{};
            HDC hdc = BeginPaint(hwnd, &paint);
            RECT rect{};
            GetClientRect(hwnd, &rect);
            FillRect(hdc, &rect, g_brush_page);
            HPEN pen = CreatePen(PS_SOLID, 1, kColBorder);
            HGDIOBJ old_pen = SelectObject(hdc, pen);
            const int center = (rect.bottom - rect.top) / 2;
            MoveToEx(hdc, rect.left + 10, center, nullptr);
            LineTo(hdc, rect.right - 10, center);
            SelectObject(hdc, old_pen);
            DeleteObject(pen);
            EndPaint(hwnd, &paint);
            return 0;
        }
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

// The tab control erases its display area with the system 3D-face colour by
// default, which reads as a grey slab behind the light child controls. Paint
// the body with the page colour so the whole tab surface stays consistent.
LRESULT CALLBACK tab_subclass_proc(
    HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam, UINT_PTR, DWORD_PTR) {
    if (message == WM_ERASEBKGND) {
        HDC hdc = reinterpret_cast<HDC>(wparam);
        RECT rect{};
        GetClientRect(hwnd, &rect);
        FillRect(hdc, &rect, g_brush_page);
        return 1;
    }
    return DefSubclassProc(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK param_panel_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_ERASEBKGND: {
            HDC hdc = reinterpret_cast<HDC>(wparam);
            RECT rect{};
            GetClientRect(hwnd, &rect);
            FillRect(hdc, &rect, g_brush_page);
            return 1;
        }
        case WM_VSCROLL: {
            RECT rect{};
            GetClientRect(hwnd, &rect);
            const int page = std::max(1, static_cast<int>(rect.bottom - rect.top));
            switch (LOWORD(wparam)) {
                case SB_LINEUP:
                    scroll_parameter_panel(-30);
                    break;
                case SB_LINEDOWN:
                    scroll_parameter_panel(30);
                    break;
                case SB_PAGEUP:
                    scroll_parameter_panel(-page);
                    break;
                case SB_PAGEDOWN:
                    scroll_parameter_panel(page);
                    break;
                case SB_THUMBTRACK:
                case SB_THUMBPOSITION: {
                    SCROLLINFO info{};
                    info.cbSize = sizeof(info);
                    info.fMask = SIF_TRACKPOS;
                    GetScrollInfo(hwnd, SB_VERT, &info);
                    scroll_parameter_panel(info.nTrackPos - g_param_scroll_y);
                    break;
                }
            }
            return 0;
        }
        case WM_MOUSEWHEEL:
            scroll_parameter_panel(-GET_WHEEL_DELTA_WPARAM(wparam) / WHEEL_DELTA * 90);
            return 0;
        case WM_SIZE:
            layout_parameter_panel();
            return 0;
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK window_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_CREATE:
            g_hwnd = hwnd;
            g_executable_dir = executable_directory();
            g_config_path = g_executable_dir / L"config.txt";
            create_fonts();
            create_controls();
            load_config();
            layout_controls();
            return 0;

        case WM_SIZE:
            layout_controls();
            return 0;

        case WM_GETMINMAXINFO: {
            auto* minmax = reinterpret_cast<MINMAXINFO*>(lparam);
            minmax->ptMinTrackSize.x = 980;
            minmax->ptMinTrackSize.y = 640;
            return 0;
        }

        // Labels, group boxes, check boxes (all "static-like" controls) and the
        // read-only log share the page surface; editable fields and combo drop
        // lists use the white field surface. This is what kills the grey wash.
        case WM_CTLCOLORSTATIC: {
            HDC dc = reinterpret_cast<HDC>(wparam);
            SetBkColor(dc, kColPage);
            SetTextColor(dc, kColText);
            return reinterpret_cast<LRESULT>(g_brush_page);
        }
        case WM_CTLCOLOREDIT:
        case WM_CTLCOLORLISTBOX: {
            HDC dc = reinterpret_cast<HDC>(wparam);
            SetBkColor(dc, kColField);
            SetTextColor(dc, kColText);
            return reinterpret_cast<LRESULT>(g_brush_field);
        }

        case WM_TIMER:
            if (wparam == kProcessTimer) {
                poll_solver();
                return 0;
            }
            break;

        case WM_NOTIFY:
            if (reinterpret_cast<NMHDR*>(lparam)->hwndFrom == g_tabs &&
                reinterpret_cast<NMHDR*>(lparam)->code == TCN_SELCHANGE) {
                show_tab(static_cast<int>(TabCtrl_GetCurSel(g_tabs)));
                layout_controls();
                return 0;
            }
            break;

        case WM_COMMAND:
            switch (LOWORD(wparam)) {
                case kSaveId:
                    save_config();
                    return 0;
                case kReloadId:
                    if (confirm_discard_changes()) load_config();
                    return 0;
                case kDefaultsId:
                    reset_defaults();
                    return 0;
                case kRunId:
                case kRunPanelId:
                    start_solver();
                    return 0;
                case kStopId:
                case kStopPanelId:
                    stop_solver();
                    return 0;
                case kOpenResultsId:
                    open_results_folder();
                    return 0;
                case kRefreshPreviewId:
                    if (save_config()) refresh_preview();
                    return 0;
                case kOpenParaviewId:
                    open_paraview_result();
                    return 0;
                case kPrevStepId:
                    move_preview_step_selection(-1);
                    return 0;
                case kNextStepId:
                    move_preview_step_selection(1);
                    return 0;
                case kSyncRawId:
                    sync_raw_from_form();
                    return 0;
                case kApplyRawId:
                    apply_raw_to_form();
                    return 0;
            }

            if (HIWORD(wparam) == EN_CHANGE ||
                HIWORD(wparam) == CBN_SELCHANGE ||
                HIWORD(wparam) == BN_CLICKED) {
                const int id = LOWORD(wparam);
                if (id >= kEditR && id <= kCheckFieldTorque &&
                    id != kEditLog &&
                    id != kEditSolverPath &&
                    id != kEditMpiRanks &&
                    id != kComboPreviewField &&
                    id != kComboPreviewStep) {
                    mark_dirty();
                }
                if ((id == kComboPreviewField || id == kComboPreviewStep) && HIWORD(wparam) == CBN_SELCHANGE && !g_loading) {
                    refresh_preview();
                }
                if (id == kComboInletType && HIWORD(wparam) == CBN_SELCHANGE) {
                    update_inlet_editor_visibility();
                }
            }
            break;

        case WM_CLOSE:
            if (g_solver.running) {
                const int result = MessageBoxW(
                    hwnd,
                    L"The solver is still running. Stop it and close?",
                    L"Pancake",
                    MB_ICONWARNING | MB_YESNO | MB_DEFBUTTON2);
                if (result != IDYES) return 0;
                stop_solver();
            }
            if (confirm_discard_changes()) DestroyWindow(hwnd);
            return 0;

        case WM_DESTROY:
            KillTimer(hwnd, kProcessTimer);
            close_solver_handles();
            if (g_tabs != nullptr) RemoveWindowSubclass(g_tabs, tab_subclass_proc, 0);
            if (g_font != nullptr) DeleteObject(g_font);
            if (g_heading_font != nullptr) DeleteObject(g_heading_font);
            if (g_mono_font != nullptr) DeleteObject(g_mono_font);
            if (g_brush_page != nullptr) DeleteObject(g_brush_page);
            if (g_brush_field != nullptr) DeleteObject(g_brush_field);
            if (g_richedit_module != nullptr) FreeLibrary(g_richedit_module);
            PostQuitMessage(0);
            return 0;
    }

    return DefWindowProcW(hwnd, message, wparam, lparam);
}
}

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show_command) {
    INITCOMMONCONTROLSEX common_controls{sizeof(INITCOMMONCONTROLSEX), ICC_TAB_CLASSES | ICC_STANDARD_CLASSES};
    InitCommonControlsEx(&common_controls);
    g_richedit_module = LoadLibraryW(L"Msftedit.dll");
    g_richedit_available = g_richedit_module != nullptr;

    g_brush_page = CreateSolidBrush(kColPage);
    g_brush_field = CreateSolidBrush(kColField);

    WNDCLASSEXW preview_class{};
    preview_class.cbSize = sizeof(preview_class);
    preview_class.lpfnWndProc = preview_proc;
    preview_class.hInstance = instance;
    preview_class.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    preview_class.hbrBackground = g_brush_field;
    preview_class.lpszClassName = kPreviewClass;
    RegisterClassExW(&preview_class);

    WNDCLASSEXW splitter_class{};
    splitter_class.cbSize = sizeof(splitter_class);
    splitter_class.lpfnWndProc = splitter_proc;
    splitter_class.hInstance = instance;
    splitter_class.hCursor = LoadCursorW(nullptr, IDC_SIZENS);
    splitter_class.hbrBackground = g_brush_page;
    splitter_class.lpszClassName = kSplitterClass;
    RegisterClassExW(&splitter_class);

    WNDCLASSEXW param_panel_class{};
    param_panel_class.cbSize = sizeof(param_panel_class);
    param_panel_class.lpfnWndProc = param_panel_proc;
    param_panel_class.hInstance = instance;
    param_panel_class.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    param_panel_class.hbrBackground = g_brush_page;
    param_panel_class.lpszClassName = kParamPanelClass;
    RegisterClassExW(&param_panel_class);

    WNDCLASSEXW window_class{};
    HICON app_icon = static_cast<HICON>(LoadImageW(
        instance, MAKEINTRESOURCEW(IDI_PANCAKE), IMAGE_ICON, 0, 0, LR_DEFAULTSIZE | LR_SHARED));
    HICON app_icon_small = static_cast<HICON>(LoadImageW(
        instance,
        MAKEINTRESOURCEW(IDI_PANCAKE),
        IMAGE_ICON,
        GetSystemMetrics(SM_CXSMICON),
        GetSystemMetrics(SM_CYSMICON),
        LR_SHARED));
    window_class.cbSize = sizeof(window_class);
    window_class.lpfnWndProc = window_proc;
    window_class.hInstance = instance;
    window_class.hIcon = app_icon != nullptr ? app_icon : LoadIconW(nullptr, IDI_APPLICATION);
    window_class.hIconSm = app_icon_small != nullptr ? app_icon_small : window_class.hIcon;
    window_class.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    window_class.hbrBackground = g_brush_page;
    window_class.lpszClassName = kWindowClass;

    if (!RegisterClassExW(&window_class)) {
        return 1;
    }

    HWND hwnd = CreateWindowExW(
        0,
        kWindowClass,
        L"Pancake",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT,
        CW_USEDEFAULT,
        1100,
        760,
        nullptr,
        nullptr,
        instance,
        nullptr);

    if (hwnd == nullptr) return 1;
    if (app_icon != nullptr) SendMessageW(hwnd, WM_SETICON, ICON_BIG, reinterpret_cast<LPARAM>(app_icon));
    if (app_icon_small != nullptr) SendMessageW(hwnd, WM_SETICON, ICON_SMALL, reinterpret_cast<LPARAM>(app_icon_small));

    ShowWindow(hwnd, show_command);
    UpdateWindow(hwnd);

    MSG message{};
    while (GetMessageW(&message, nullptr, 0, 0) > 0) {
        TranslateMessage(&message);
        DispatchMessageW(&message);
    }

    return static_cast<int>(message.wParam);
}
