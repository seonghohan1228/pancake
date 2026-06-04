#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <windows.h>
#include <commctrl.h>
#include <commdlg.h>
#include <richedit.h>
#include <shellapi.h>
#include <windowsx.h>

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
#include <set>
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
constexpr wchar_t kRailClass[] = L"PancakeRailPanel";
constexpr wchar_t kSummaryClass[] = L"PancakeSummaryCard";
constexpr UINT_PTR kProcessTimer = 1;
constexpr double kPi = 3.14159265358979323846264338327950288;

// Layout metrics (device-independent-ish; window has a sensible minimum size).
constexpr int kActionBarH = 56;
constexpr int kMargin = 10;
constexpr int kGap = 8;
constexpr int kSplitterW = 6;
constexpr int kRailToolbarH = 40;
constexpr int kSummaryH = 156;
constexpr int kRailMinWidth = 300;
constexpr int kRailMaxSlack = 460;   // right area must keep at least this many px
constexpr int kConsoleMinHeight = 90;
constexpr int kViewportMinHeight = 220;

// Modern light palette. The whole surface shares one page colour so labels,
// section headers and the summary card read as a single clean canvas instead of
// the default system grey.
constexpr COLORREF kColPage = RGB(244, 247, 251);   // window / panel surface
constexpr COLORREF kColCard = RGB(255, 255, 255);   // raised cards / canvas
constexpr COLORREF kColField = RGB(255, 255, 255);  // edit & list backgrounds
constexpr COLORREF kColText = RGB(28, 36, 46);      // primary text
constexpr COLORREF kColMuted = RGB(108, 119, 133);  // secondary text
constexpr COLORREF kColBorder = RGB(208, 215, 224); // hairline separators
constexpr COLORREF kColAccent = RGB(38, 110, 196);  // headings / highlights
constexpr COLORREF kColRun = RGB(33, 140, 84);      // run button
constexpr COLORREF kColRunHot = RGB(40, 160, 98);
constexpr COLORREF kColStop = RGB(176, 58, 63);     // stop button
constexpr COLORREF kColStopHot = RGB(196, 70, 75);
constexpr COLORREF kColInvalidBg = RGB(255, 236, 236);
constexpr COLORREF kColInvalidText = RGB(176, 40, 40);
constexpr COLORREF kColWarn = RGB(150, 100, 15);
constexpr COLORREF kColOk = RGB(33, 140, 84);

enum ControlId {
    // Action bar / files / run
    kRunId = 1000,
    kStopId,
    kSaveId,
    kSaveAsId,
    kOpenId,
    kReloadId,
    kDefaultsId,
    kSolverBrowseId,
    kComboPreset,
    kComboUnitSystem,
    kProgressId,
    kEditMpiRanks,

    // Rail toolbar
    kEditSearch,
    kCheckShowAdvanced,

    // Results toolbar / console
    kComboPreviewField,
    kComboPreviewStep,
    kPrevStepId,
    kNextStepId,
    kRefreshPreviewId,
    kOpenResultsId,
    kOpenParaviewId,
    kClearLogId,
    kPreviewCanvasId,
    kEditLog,

    // Raw config
    kSyncRawId,
    kApplyRawId,
    kEditRawConfig,

    // Parameter edits
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
    kEditInletTheta,
    kEditInletZ,
    kEditInletSize,
    kEditInletPressure,
    kEditOutputDir,
    kEditFilenamePrefix,
    kEditExternalLoadMagnitude,
    kEditExternalLoadDirection,
    kEditExternalLoadZ,
    kEditBearingInitialX,
    kEditBearingInitialY,
    kEditBearingInitialZ,
    kEditBearingInitialVx,
    kEditBearingInitialVy,
    kEditBearingInitialVz,
    kEditBearingMass,
    kEditBearingStiffnessX,
    kEditBearingStiffnessY,
    kEditBearingStiffnessZ,
    kEditBearingDampingX,
    kEditBearingDampingY,
    kEditBearingDampingZ,
    kEditMinFilmThickness,
    kEditTemperatureInitial,
    kEditTemperatureReference,
    kEditJournalWallTemperature,
    kEditBearingWallTemperature,
    kEditRhoCp,
    kEditThermalConductivity,
    kEditJournalHeatTransfer,
    kEditBearingHeatTransfer,
    kEditLastParameter = kEditBearingHeatTransfer,

    // Combos
    kComboCavitation = 3000,
    kComboBcSouth,
    kComboBcNorth,
    kComboInletType,
    kComboSolutionMode,
    kComboMotionModel,
    kComboTemperatureModel,
    kComboPressureTimeMethod,
    kComboMotionTimeMethod,
    kComboTemperatureTimeMethod,
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
    kComboUnitPCav,
    kComboUnitBulkModulus,
    kComboUnitBcSouthVal,
    kComboUnitBcNorthVal,
    kComboUnitInletTheta,
    kComboUnitInletZ,
    kComboUnitInletSize,
    kComboUnitInletPressure,

    // Output field checkboxes
    kCheckLogOuter = 4000,
    kCheckBearingInitialFromAttitude,
    kCheckStopOnNonpositiveFilm,
    kCheckOutput3d,
    kCheckOutputFlat,
    kCheckFieldPressure,
    kCheckFieldFilmContent,
    kCheckFieldH,
    kCheckFieldRho,
    kCheckFieldTemperature,
    kCheckFieldHeatGeneration,
    kCheckFieldInlet,
    kCheckFieldVelocity,
    kCheckFieldPressureForceX,
    kCheckFieldPressureForceY,
    kCheckFieldPressureForceZ,
    kCheckFieldLoadX,
    kCheckFieldLoadY,
    kCheckFieldLoadZ,
    kCheckFieldViscousX,
    kCheckFieldViscousY,
    kCheckFieldViscousZ,
    kCheckFieldFluidX,
    kCheckFieldFluidY,
    kCheckFieldFluidZ,
    kCheckFieldExternalX,
    kCheckFieldExternalY,
    kCheckFieldExternalZ,
    kCheckFieldBearingX,
    kCheckFieldBearingY,
    kCheckFieldBearingZ,
    kCheckFieldBearingAttitude,
    kCheckFieldTorque
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

enum class ParamVisibility { Always, CircularInletOnly };
enum class RowKind { Edit, Combo, Check, RawEditor };

// One input row in the rail. Labels and fixed-unit suffixes are painted by the
// rail panel (no STATIC child controls), which keeps the child-control count low
// and removes the combo/group-box z-order fights that plagued the old GUI.
struct Row {
    RowKind kind = RowKind::Edit;
    int section = 0;
    bool advanced = false;
    ParamVisibility visibility = ParamVisibility::Always;
    int value_id = 0;       // edit / combo / checkbox id
    int unit_id = 0;        // override unit combo id, or 0 for a painted suffix
    std::wstring label;     // painted label (left column)
    std::wstring suffix;    // painted unit suffix when unit_id == 0
    HWND value = nullptr;
    HWND unit = nullptr;
    // runtime layout
    bool visible = false;
    RECT label_rect{};
    RECT suffix_rect{};
};

struct Section {
    std::wstring title;
    bool advanced = false;  // only the header's "all rows advanced" hint
    bool collapsed = false;
    bool is_raw = false;
    // runtime
    bool visible = false;
    RECT header_rect{};
};

struct PreviewData {
    bool ok = false;
    int nx = 0;
    int ny = 0;
    int step = -1;
    double min_value = 0.0;
    double max_value = 0.0;
    double mean_value = 0.0;
    std::string field;
    std::wstring message = L"No result loaded.";
    std::wstring source_path;
    std::vector<double> values;
};

struct PreviewStep {
    int step = -1;
    fs::path root;
    fs::file_time_type write_time{};
    std::vector<fs::path> files;
};

struct SolverProcess {
    PROCESS_INFORMATION process_info{};
    HANDLE stdout_read = nullptr;
    bool running = false;
};

enum class RunState { Idle, Running, Done, Failed };

const std::array<OutputFieldControl, 28> kOutputFields = {{
    {kCheckFieldPressure, "pressure", L"pressure"},
    {kCheckFieldFilmContent, "film_content", L"film_content"},
    {kCheckFieldH, "h", L"h (film thickness)"},
    {kCheckFieldRho, "rho", L"rho"},
    {kCheckFieldTemperature, "temperature", L"temperature"},
    {kCheckFieldHeatGeneration, "heat_generation", L"heat_generation"},
    {kCheckFieldInlet, "inlet_indicator", L"inlet_indicator"},
    {kCheckFieldVelocity, "velocity", L"velocity"},
    {kCheckFieldPressureForceX, "pressure_force_x", L"pressure_force_x"},
    {kCheckFieldPressureForceY, "pressure_force_y", L"pressure_force_y"},
    {kCheckFieldPressureForceZ, "pressure_force_z", L"pressure_force_z"},
    {kCheckFieldLoadX, "load_x", L"load_x (legacy pressure)"},
    {kCheckFieldLoadY, "load_y", L"load_y (legacy pressure)"},
    {kCheckFieldLoadZ, "load_z", L"load_z (legacy pressure)"},
    {kCheckFieldViscousX, "viscous_force_x", L"viscous_force_x"},
    {kCheckFieldViscousY, "viscous_force_y", L"viscous_force_y"},
    {kCheckFieldViscousZ, "viscous_force_z", L"viscous_force_z"},
    {kCheckFieldFluidX, "fluid_force_x", L"fluid_force_x"},
    {kCheckFieldFluidY, "fluid_force_y", L"fluid_force_y"},
    {kCheckFieldFluidZ, "fluid_force_z", L"fluid_force_z"},
    {kCheckFieldExternalX, "external_load_x", L"external_load_x"},
    {kCheckFieldExternalY, "external_load_y", L"external_load_y"},
    {kCheckFieldExternalZ, "external_load_z", L"external_load_z"},
    {kCheckFieldBearingX, "bearing_x", L"bearing_x"},
    {kCheckFieldBearingY, "bearing_y", L"bearing_y"},
    {kCheckFieldBearingZ, "bearing_z", L"bearing_z"},
    {kCheckFieldBearingAttitude, "bearing_attitude_angle", L"bearing_attitude_angle"},
    {kCheckFieldTorque, "friction_torque", L"friction_torque"},
}};

// ----- globals ------------------------------------------------------------
HWND g_hwnd = nullptr;
HWND g_rail = nullptr;          // scrollable input panel
HWND g_summary = nullptr;       // pinned summary/validation card
HWND g_preview_canvas = nullptr;
HWND g_rail_splitter = nullptr; // vertical, resizes rail width
HWND g_console_splitter = nullptr; // horizontal, resizes console height
HWND g_run_button = nullptr;
HWND g_stop_button = nullptr;
HWND g_progress = nullptr;
HWND g_parent_override = nullptr;

HFONT g_font = nullptr;
HFONT g_heading_font = nullptr;
HFONT g_small_font = nullptr;
HFONT g_mono_font = nullptr;
HBRUSH g_brush_page = nullptr;
HBRUSH g_brush_field = nullptr;
HBRUSH g_brush_card = nullptr;
HBRUSH g_brush_invalid = nullptr;
HMODULE g_richedit_module = nullptr;
bool g_richedit_available = false;

fs::path g_executable_dir;
fs::path g_config_path;
std::wstring g_solver_path;     // optional explicit solver exe override
SimulationConfig g_config;
std::map<std::string, std::string> g_config_value_text;
std::string g_config_inlets_text;
SolverProcess g_solver;
PreviewData g_preview;
std::vector<PreviewStep> g_preview_steps;

bool g_dirty = false;
bool g_loading = false;
int g_rail_width = 360;
int g_console_height = 150;
int g_rail_scroll_y = 0;
int g_rail_content_height = 0;
bool g_rail_splitter_dragging = false;
bool g_console_splitter_dragging = false;

RunState g_run_state = RunState::Idle;
ULONGLONG g_run_start_tick = 0;
double g_run_progress = 0.0;     // 0..1 parsed from solver stdout
double g_run_elapsed_s = 0.0;

std::vector<ControlRef> g_control_refs;
std::vector<Section> g_sections;
std::vector<Row> g_rows;
std::set<int> g_invalid_ids;
std::vector<std::pair<std::wstring, std::wstring>> g_summary_pairs; // label,value
std::vector<std::wstring> g_problems;   // validation messages (errors + warnings)
bool g_has_errors = false;

// ----- string / path / file helpers --------------------------------------
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
    std::wstring out;
    out.reserve(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] == L'\r') {
            out += L"\r\n";
            if (i + 1 < text.size() && text[i + 1] == L'\n') ++i;
        } else if (text[i] == L'\n') {
            out += L"\r\n";
        } else {
            out += text[i];
        }
    }
    return out;
}

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

// ----- control registry / creation ---------------------------------------
HWND find_control(int id) {
    for (const auto& ref : g_control_refs) {
        if (ref.id == id) return ref.hwnd;
    }
    return nullptr;
}

void remember_control(int id, HWND hwnd) {
    if (id != 0) g_control_refs.push_back({id, hwnd});
}

void move_control(HWND hwnd, int x, int y, int width, int height, BOOL repaint = TRUE) {
    if (hwnd != nullptr) MoveWindow(hwnd, x, y, std::max(width, 1), std::max(height, 1), repaint);
}

void move_control_id(int id, int x, int y, int width, int height, BOOL repaint = TRUE) {
    move_control(find_control(id), x, y, width, height, repaint);
}

HWND create_child(HWND parent, const wchar_t* cls, const wchar_t* text,
                  DWORD style, DWORD ex_style, int id, HFONT font = nullptr) {
    HWND hwnd = CreateWindowExW(
        ex_style, cls, text, WS_CHILD | WS_CLIPSIBLINGS | style,
        0, 0, 10, 10, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr), nullptr);
    SendMessageW(hwnd, WM_SETFONT, reinterpret_cast<WPARAM>(font ? font : g_font), TRUE);
    remember_control(id, hwnd);
    return hwnd;
}

HWND create_edit(HWND parent, int id, DWORD extra = ES_AUTOHSCROLL) {
    return create_child(parent, L"EDIT", L"", WS_TABSTOP | extra, WS_EX_CLIENTEDGE, id);
}

HWND create_button(HWND parent, int id, const wchar_t* text, DWORD extra = 0) {
    return create_child(parent, L"BUTTON", text,
                        WS_VISIBLE | WS_TABSTOP | BS_PUSHBUTTON | extra, 0, id);
}

HWND create_checkbox(HWND parent, int id, const wchar_t* text) {
    return create_child(parent, L"BUTTON", text, WS_TABSTOP | BS_AUTOCHECKBOX, 0, id);
}

HWND create_combo(HWND parent, int id, const std::vector<const wchar_t*>& items) {
    HWND combo = CreateWindowExW(
        0, L"COMBOBOX", L"",
        WS_CHILD | WS_CLIPSIBLINGS | WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST | CBS_HASSTRINGS,
        0, 0, 10, 240, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr), nullptr);
    SendMessageW(combo, WM_SETFONT, reinterpret_cast<WPARAM>(g_font), TRUE);
    for (const wchar_t* item : items) {
        SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(item));
    }
    if (!items.empty()) SendMessageW(combo, CB_SETCURSEL, 0, 0);
    remember_control(id, combo);
    return combo;
}

// ----- small accessors -----------------------------------------------------
void set_edit_text(int id, const std::string& text) {
    HWND hwnd = find_control(id);
    if (hwnd != nullptr) SetWindowTextW(hwnd, normalize_newlines_for_edit(utf8_to_wide(text)).c_str());
}
void set_edit_text(int id, double value) { set_edit_text(id, compact_number(value)); }
void set_edit_text(int id, int value) { set_edit_text(id, std::to_string(value)); }

std::string get_edit_utf8(int id) {
    HWND hwnd = find_control(id);
    if (hwnd == nullptr) return {};
    return wide_to_utf8(get_window_text(hwnd));
}

bool combo_has_selection(int id) {
    HWND hwnd = find_control(id);
    return hwnd != nullptr && SendMessageW(hwnd, CB_GETCURSEL, 0, 0) >= 0;
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
    if (!combo_has_selection(id) && !items.empty()) SendMessageW(hwnd, CB_SETCURSEL, 0, 0);
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

bool combo_is(int id, const wchar_t* value) { return get_combo_text(id) == value; }

void mark_dirty();
void recompute_summary();
void layout_everything();
void layout_rail();
void update_title();

// ----- unit conversion ----------------------------------------------------
double unit_scale_to_meters(int unit_id) {
    const std::wstring unit = get_combo_text(unit_id);
    if (unit == L"cm") return 1.0e-2;
    if (unit == L"mm") return 1.0e-3;
    if (unit == L"um") return 1.0e-6;
    return 1.0;
}
double unit_scale_to_seconds(int unit_id) {
    return get_combo_text(unit_id) == L"ms" ? 1.0e-3 : 1.0;
}
double unit_scale_to_pascals(int unit_id) {
    const std::wstring unit = get_combo_text(unit_id);
    if (unit == L"kPa") return 1.0e3;
    if (unit == L"MPa") return 1.0e6;
    return 1.0;
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

// Validating reads (used on Save/Run/Apply): they raise a MessageBox on bad
// input so the user gets a precise complaint.
bool read_double(int id, const wchar_t* label, double& value) {
    const std::string text = trim_copy(get_edit_utf8(id));
    try {
        size_t used = 0;
        value = std::stod(text, &used);
        if (used != text.size()) throw std::invalid_argument("trailing");
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
        if (used != text.size()) throw std::invalid_argument("trailing");
        return true;
    } catch (...) {
        MessageBoxW(g_hwnd, (std::wstring(L"Invalid integer for ") + label).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
}

// Silent reads (used by the live summary/validation): no popups.
bool try_double(int id, double& value) {
    const std::string text = trim_copy(get_edit_utf8(id));
    if (text.empty()) return false;
    try {
        size_t used = 0;
        value = std::stod(text, &used);
        return used == text.size();
    } catch (...) { return false; }
}
bool try_int(int id, int& value) {
    const std::string text = trim_copy(get_edit_utf8(id));
    if (text.empty()) return false;
    try {
        size_t used = 0;
        value = std::stoi(text, &used);
        return used == text.size();
    } catch (...) { return false; }
}

double si_length(int id, int unit_id) {
    double v = 0.0;
    if (!try_double(id, v)) return std::numeric_limits<double>::quiet_NaN();
    return v * unit_scale_to_meters(unit_id);
}
double si_omega() {
    double v = 0.0;
    if (!try_double(kEditOmega, v)) return std::numeric_limits<double>::quiet_NaN();
    return combo_is(kComboUnitOmega, L"rpm") ? v * 2.0 * kPi / 60.0 : v;
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
    if (combo_is(unit_id, L"deg")) slope = std::tan(value * kPi / 180.0);
    else if (combo_is(unit_id, L"rad")) slope = std::tan(value);
    else slope = value;
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
CavitationModel cavitation_from_combo() {
    const std::wstring value = get_combo_text(kComboCavitation);
    if (value == L"GUMBEL") return CavitationModel::GUMBEL;
    if (value == L"FULL_SOMMERFELD") return CavitationModel::FULL_SOMMERFELD;
    return CavitationModel::ELROD_ADAMS;
}

void set_motion_model_combo(MotionModel model) {
    set_combo_text(kComboMotionModel, utf8_to_wide(to_config_value(model)).c_str());
}
MotionModel motion_model_from_combo() {
    return motion_model_from_config(wide_to_utf8(get_combo_text(kComboMotionModel)));
}
void set_time_method_combo(int id, TimeSteppingMethod method) {
    set_combo_text(id, utf8_to_wide(to_config_value(method)).c_str());
}
TimeSteppingMethod time_method_from_combo(int id) {
    return time_stepping_from_config(wide_to_utf8(get_combo_text(id)));
}
void set_solution_mode_combo(SolutionMode mode) {
    set_combo_text(kComboSolutionMode, utf8_to_wide(to_config_value(mode)).c_str());
}
SolutionMode solution_mode_from_combo() {
    return solution_mode_from_config(wide_to_utf8(get_combo_text(kComboSolutionMode)));
}
void set_temperature_model_combo(TemperatureModel model) {
    set_combo_text(kComboTemperatureModel, utf8_to_wide(to_config_value(model)).c_str());
}
TemperatureModel temperature_model_from_combo() {
    return temperature_model_from_config(wide_to_utf8(get_combo_text(kComboTemperatureModel)));
}

std::string native_or_converted_text(int edit_id, int unit_id, double native_value, bool native_unit) {
    return native_unit ? trim_copy(get_edit_utf8(edit_id)) : compact_number(native_value);
}

// ----- config <-> form ----------------------------------------------------
void set_edit_from_config_text(int id, const char* key, double fallback) {
    auto it = g_config_value_text.find(key);
    set_edit_text(id, it != g_config_value_text.end() ? it->second : compact_number(fallback));
}
void set_edit_from_config_text(int id, const char* key, int fallback) {
    auto it = g_config_value_text.find(key);
    set_edit_text(id, it != g_config_value_text.end() ? it->second : std::to_string(fallback));
}
void set_length_from_config_text(int edit_id, int unit_id, const char* key, double meters) {
    if (unit_is_native_length(unit_id)) set_edit_from_config_text(edit_id, key, meters);
    else set_edit_text(edit_id, compact_number(meters / unit_scale_to_meters(unit_id)));
}
void set_time_from_config_text(int edit_id, int unit_id, const char* key, double seconds) {
    if (unit_is_native_time(unit_id)) set_edit_from_config_text(edit_id, key, seconds);
    else set_edit_text(edit_id, compact_number(seconds / unit_scale_to_seconds(unit_id)));
}
void set_pressure_from_config_text(int edit_id, int unit_id, const char* key, double pascals) {
    if (unit_is_native_pressure(unit_id)) set_edit_from_config_text(edit_id, key, pascals);
    else set_edit_text(edit_id, compact_number(pascals / unit_scale_to_pascals(unit_id)));
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
    if (!combo_has_selection(kComboUnitPCav)) set_combo_text(kComboUnitPCav, L"Pa");
    if (!combo_has_selection(kComboUnitBulkModulus)) set_combo_text(kComboUnitBulkModulus, L"Pa");
    if (!combo_has_selection(kComboUnitBcSouthVal)) set_combo_text(kComboUnitBcSouthVal, L"Pa");
    if (!combo_has_selection(kComboUnitBcNorthVal)) set_combo_text(kComboUnitBcNorthVal, L"Pa");
    if (!combo_has_selection(kComboUnitInletTheta)) set_combo_text(kComboUnitInletTheta, L"deg");
    if (!combo_has_selection(kComboUnitInletZ)) set_combo_text(kComboUnitInletZ, L"m");
    if (!combo_has_selection(kComboUnitInletSize)) set_combo_text(kComboUnitInletSize, L"deg");
    if (!combo_has_selection(kComboUnitInletPressure)) set_combo_text(kComboUnitInletPressure, L"Pa");
}

struct InletTextFields {
    bool valid = false;
    InletConfig::Type type = InletConfig::Type::GROOVE;
    std::string theta, z, size, p_supply;
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
    if (combo_is(kComboUnitInletTheta, L"rad"))
        set_edit_text(kEditInletTheta, compact_number(inlet.theta * kPi / 180.0));
    else
        set_edit_text(kEditInletTheta, use_raw ? raw.theta : compact_number(inlet.theta));
    set_edit_text(kEditInletZ, (use_raw && unit_is_native_length(kComboUnitInletZ)) ? raw.z : compact_number(inlet.z / unit_scale_to_meters(kComboUnitInletZ)));
    if (inlet.type == InletConfig::Type::CIRCULAR) {
        set_edit_text(kEditInletSize, (use_raw && unit_is_native_length(kComboUnitInletSize)) ? raw.size : compact_number(inlet.size / unit_scale_to_meters(kComboUnitInletSize)));
    } else if (combo_is(kComboUnitInletSize, L"rad")) {
        set_edit_text(kEditInletSize, compact_number(inlet.size * kPi / 180.0));
    } else {
        set_edit_text(kEditInletSize, use_raw ? raw.size : compact_number(inlet.size));
    }
    set_edit_text(kEditInletPressure,
        (use_raw && unit_is_native_pressure(kComboUnitInletPressure)) ? raw.p_supply
        : compact_number(inlet.p_supply / unit_scale_to_pascals(kComboUnitInletPressure)));
}

// The inlet z row (circular only) and the size unit set + label switch with
// the inlet type. The row's runtime visibility is handled by the rail layout
// via ParamVisibility::CircularInletOnly; here we only update the combo set,
// the painted label text, and trigger a relayout.
std::wstring g_inlet_size_label = L"width";
void update_inlet_editor_visibility() {
    const bool circular = combo_is(kComboInletType, L"CIRCULAR");
    set_combo_items(kComboUnitInletSize,
        circular ? std::vector<const wchar_t*>{L"m", L"mm", L"cm", L"um"}
                 : std::vector<const wchar_t*>{L"deg", L"rad"},
        circular ? L"mm" : L"deg");
    g_inlet_size_label = circular ? L"radius" : L"width";
    for (Row& row : g_rows) {
        if (row.value_id == kEditInletSize) row.label = g_inlet_size_label;
    }
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

std::string form_config_text(const SimulationConfig& cfg);

void populate_form_from_config(const SimulationConfig& cfg) {
    g_loading = true;
    set_native_unit_defaults();

    set_length_from_config_text(kEditR, kComboUnitR, "R", cfg.R);
    set_length_from_config_text(kEditC, kComboUnitC, "c", cfg.c);
    set_length_from_config_text(kEditE, kComboUnitE, "e", cfg.e);
    set_length_from_config_text(kEditL, kComboUnitL, "L", cfg.L);
    if (combo_is(kComboUnitAttitude, L"rad")) set_edit_text(kEditAttitude, compact_number(cfg.attitude_angle_deg * kPi / 180.0));
    else set_edit_from_config_text(kEditAttitude, "attitude_angle_deg", cfg.attitude_angle_deg);
    if (combo_is(kComboUnitLoadAngle, L"rad")) set_edit_text(kEditLoadAngle, compact_number(cfg.load_angle_deg * kPi / 180.0));
    else set_edit_from_config_text(kEditLoadAngle, "load_angle_deg", cfg.load_angle_deg);
    if (combo_is(kComboUnitTiltX, L"deg")) set_edit_text(kEditTiltX, compact_number(std::atan(cfg.tilt_slope_x) * 180.0 / kPi));
    else if (combo_is(kComboUnitTiltX, L"rad")) set_edit_text(kEditTiltX, compact_number(std::atan(cfg.tilt_slope_x)));
    else set_edit_from_config_text(kEditTiltX, "tilt_slope_x", cfg.tilt_slope_x);
    if (combo_is(kComboUnitTiltY, L"deg")) set_edit_text(kEditTiltY, compact_number(std::atan(cfg.tilt_slope_y) * 180.0 / kPi));
    else if (combo_is(kComboUnitTiltY, L"rad")) set_edit_text(kEditTiltY, compact_number(std::atan(cfg.tilt_slope_y)));
    else set_edit_from_config_text(kEditTiltY, "tilt_slope_y", cfg.tilt_slope_y);

    set_edit_from_config_text(kEditNTheta, "n_theta_global", cfg.n_theta_global);
    set_edit_from_config_text(kEditNZ, "n_z_global", cfg.n_z_global);
    set_time_from_config_text(kEditEndT, kComboUnitEndT, "end_t", cfg.end_t);
    set_time_from_config_text(kEditDt, kComboUnitDt, "dt", cfg.dt);
    set_time_from_config_text(kEditWriteInterval, kComboUnitWriteInterval, "write_interval", cfg.write_interval);
    set_solution_mode_combo(cfg.solution_mode);
    set_time_method_combo(kComboPressureTimeMethod, cfg.pressure_time_method);
    set_time_method_combo(kComboMotionTimeMethod, cfg.motion_time_method);
    set_time_method_combo(kComboTemperatureTimeMethod, cfg.temperature_time_method);
    if (combo_is(kComboUnitOmega, L"rpm")) set_edit_text(kEditOmega, compact_number(cfg.omega * 60.0 / (2.0 * kPi)));
    else set_edit_from_config_text(kEditOmega, "omega", cfg.omega);
    set_edit_from_config_text(kEditMu, "mu", cfg.mu);
    set_edit_from_config_text(kEditRho, "rho", cfg.rho);
    set_pressure_from_config_text(kEditPCav, kComboUnitPCav, "p_cav", cfg.p_cav);
    set_pressure_from_config_text(kEditBulkModulus, kComboUnitBulkModulus, "bulk_modulus", cfg.bulk_modulus);
    set_temperature_model_combo(cfg.temperature_model);
    set_edit_from_config_text(kEditTemperatureInitial, "temperature_initial", cfg.temperature_initial);
    set_edit_from_config_text(kEditTemperatureReference, "temperature_reference", cfg.temperature_reference);
    set_edit_from_config_text(kEditJournalWallTemperature, "journal_wall_temperature", cfg.journal_wall_temperature);
    set_edit_from_config_text(kEditBearingWallTemperature, "bearing_wall_temperature", cfg.bearing_wall_temperature);
    set_edit_from_config_text(kEditRhoCp, "rho_cp", cfg.rho_cp);
    set_edit_from_config_text(kEditThermalConductivity, "thermal_conductivity", cfg.thermal_conductivity);
    set_edit_from_config_text(kEditJournalHeatTransfer, "journal_heat_transfer", cfg.journal_heat_transfer);
    set_edit_from_config_text(kEditBearingHeatTransfer, "bearing_heat_transfer", cfg.bearing_heat_transfer);
    set_motion_model_combo(cfg.motion_model);
    const double external_load_magnitude = std::hypot(cfg.external_load_x, cfg.external_load_y);
    const double external_load_direction = external_load_magnitude > 0.0
        ? std::atan2(cfg.external_load_y, cfg.external_load_x) * 180.0 / kPi
        : 0.0;
    set_edit_text(kEditExternalLoadMagnitude, external_load_magnitude);
    set_edit_text(kEditExternalLoadDirection, external_load_direction);
    set_edit_from_config_text(kEditExternalLoadZ, "external_load_z", cfg.external_load_z);
    set_check(kCheckBearingInitialFromAttitude, cfg.bearing_initial_from_attitude);
    set_edit_from_config_text(kEditBearingInitialX, "bearing_initial_x", cfg.bearing_initial_x);
    set_edit_from_config_text(kEditBearingInitialY, "bearing_initial_y", cfg.bearing_initial_y);
    set_edit_from_config_text(kEditBearingInitialZ, "bearing_initial_z", cfg.bearing_initial_z);
    set_edit_from_config_text(kEditBearingInitialVx, "bearing_initial_vx", cfg.bearing_initial_vx);
    set_edit_from_config_text(kEditBearingInitialVy, "bearing_initial_vy", cfg.bearing_initial_vy);
    set_edit_from_config_text(kEditBearingInitialVz, "bearing_initial_vz", cfg.bearing_initial_vz);
    set_edit_from_config_text(kEditBearingMass, "bearing_mass", cfg.bearing_mass);
    set_edit_from_config_text(kEditBearingStiffnessX, "bearing_stiffness_x", cfg.bearing_stiffness_x);
    set_edit_from_config_text(kEditBearingStiffnessY, "bearing_stiffness_y", cfg.bearing_stiffness_y);
    set_edit_from_config_text(kEditBearingStiffnessZ, "bearing_stiffness_z", cfg.bearing_stiffness_z);
    set_edit_from_config_text(kEditBearingDampingX, "bearing_damping_x", cfg.bearing_damping_x);
    set_edit_from_config_text(kEditBearingDampingY, "bearing_damping_y", cfg.bearing_damping_y);
    set_edit_from_config_text(kEditBearingDampingZ, "bearing_damping_z", cfg.bearing_damping_z);
    set_edit_from_config_text(kEditMinFilmThickness, "min_film_thickness", cfg.min_film_thickness);
    set_check(kCheckStopOnNonpositiveFilm, cfg.stop_on_nonpositive_film);
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
    for (const auto& field : kOutputFields) set_check(field.id, cfg.output_field_enabled(field.key));

    if (g_solver_path.empty()) g_solver_path = (g_executable_dir / L"pancake.exe").wstring();
    if (trim_copy(get_edit_utf8(kEditMpiRanks)).empty()) set_edit_text(kEditMpiRanks, 1);
    set_edit_text(kEditRawConfig, form_config_text(cfg));
    g_loading = false;
    if (g_rail != nullptr) layout_rail();
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
    cfg.solution_mode = solution_mode_from_combo();
    cfg.pressure_time_method = time_method_from_combo(kComboPressureTimeMethod);
    cfg.motion_time_method = time_method_from_combo(kComboMotionTimeMethod);
    cfg.temperature_time_method = time_method_from_combo(kComboTemperatureTimeMethod);
    if (!read_omega_rad_s(kEditOmega, kComboUnitOmega, L"omega", cfg.omega)) return false;
    if (!read_double(kEditMu, L"mu", cfg.mu)) return false;
    if (!read_double(kEditRho, L"rho", cfg.rho)) return false;
    if (!read_pressure_pascals(kEditPCav, kComboUnitPCav, L"p_cav", cfg.p_cav)) return false;
    if (!read_pressure_pascals(kEditBulkModulus, kComboUnitBulkModulus, L"bulk_modulus", cfg.bulk_modulus)) return false;
    cfg.temperature_model = temperature_model_from_combo();
    if (!read_double(kEditTemperatureInitial, L"temperature_initial", cfg.temperature_initial)) return false;
    if (!read_double(kEditTemperatureReference, L"temperature_reference", cfg.temperature_reference)) return false;
    if (!read_double(kEditJournalWallTemperature, L"journal_wall_temperature", cfg.journal_wall_temperature)) return false;
    if (!read_double(kEditBearingWallTemperature, L"bearing_wall_temperature", cfg.bearing_wall_temperature)) return false;
    if (!read_double(kEditRhoCp, L"rho_cp", cfg.rho_cp)) return false;
    if (!read_double(kEditThermalConductivity, L"thermal_conductivity", cfg.thermal_conductivity)) return false;
    if (!read_double(kEditJournalHeatTransfer, L"journal_heat_transfer", cfg.journal_heat_transfer)) return false;
    if (!read_double(kEditBearingHeatTransfer, L"bearing_heat_transfer", cfg.bearing_heat_transfer)) return false;
    cfg.motion_model = motion_model_from_combo();
    double external_load_magnitude = 0.0;
    double external_load_direction_deg = 0.0;
    if (!read_double(kEditExternalLoadMagnitude, L"external_load_magnitude", external_load_magnitude)) return false;
    if (!read_double(kEditExternalLoadDirection, L"external_load_direction", external_load_direction_deg)) return false;
    const double external_load_direction_rad = external_load_direction_deg * kPi / 180.0;
    cfg.external_load_x = external_load_magnitude * std::cos(external_load_direction_rad);
    cfg.external_load_y = external_load_magnitude * std::sin(external_load_direction_rad);
    if (!read_double(kEditExternalLoadZ, L"external_load_z", cfg.external_load_z)) return false;
    cfg.bearing_initial_from_attitude = get_check(kCheckBearingInitialFromAttitude);
    if (!read_double(kEditBearingInitialX, L"bearing_initial_x", cfg.bearing_initial_x)) return false;
    if (!read_double(kEditBearingInitialY, L"bearing_initial_y", cfg.bearing_initial_y)) return false;
    if (!read_double(kEditBearingInitialZ, L"bearing_initial_z", cfg.bearing_initial_z)) return false;
    if (!read_double(kEditBearingInitialVx, L"bearing_initial_vx", cfg.bearing_initial_vx)) return false;
    if (!read_double(kEditBearingInitialVy, L"bearing_initial_vy", cfg.bearing_initial_vy)) return false;
    if (!read_double(kEditBearingInitialVz, L"bearing_initial_vz", cfg.bearing_initial_vz)) return false;
    if (!read_double(kEditBearingMass, L"bearing_mass", cfg.bearing_mass)) return false;
    if (!read_double(kEditBearingStiffnessX, L"bearing_stiffness_x", cfg.bearing_stiffness_x)) return false;
    if (!read_double(kEditBearingStiffnessY, L"bearing_stiffness_y", cfg.bearing_stiffness_y)) return false;
    if (!read_double(kEditBearingStiffnessZ, L"bearing_stiffness_z", cfg.bearing_stiffness_z)) return false;
    if (!read_double(kEditBearingDampingX, L"bearing_damping_x", cfg.bearing_damping_x)) return false;
    if (!read_double(kEditBearingDampingY, L"bearing_damping_y", cfg.bearing_damping_y)) return false;
    if (!read_double(kEditBearingDampingZ, L"bearing_damping_z", cfg.bearing_damping_z)) return false;
    if (!read_double(kEditMinFilmThickness, L"min_film_thickness", cfg.min_film_thickness)) return false;
    cfg.stop_on_nonpositive_film = get_check(kCheckStopOnNonpositiveFilm);
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
    for (const auto& field : kOutputFields) if (get_check(field.id)) cfg.output_fields.push_back(field.key);

    if (cfg.output_dir.empty()) { MessageBoxW(g_hwnd, L"output_dir cannot be empty.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    if (cfg.filename_prefix.empty()) { MessageBoxW(g_hwnd, L"filename_prefix cannot be empty.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    if (cfg.n_theta_global <= 0 || cfg.n_z_global <= 0) { MessageBoxW(g_hwnd, L"Grid sizes must be positive.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    if (cfg.dt <= 0.0 || cfg.write_interval <= 0.0) { MessageBoxW(g_hwnd, L"dt and write_interval must be positive.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    if (cfg.temperature_model == TemperatureModel::ENERGY_EQUATION) {
        if (cfg.rho_cp <= 0.0 || cfg.thermal_conductivity <= 0.0) {
            MessageBoxW(g_hwnd, L"rho_cp and thermal_conductivity must be positive for ENERGY_EQUATION.", L"Pancake", MB_ICONERROR | MB_OK);
            return false;
        }
        if (cfg.journal_heat_transfer < 0.0 || cfg.bearing_heat_transfer < 0.0) {
            MessageBoxW(g_hwnd, L"Thermal wall heat-transfer coefficients cannot be negative.", L"Pancake", MB_ICONERROR | MB_OK);
            return false;
        }
    }
    if (!std::isfinite(external_load_magnitude) || !std::isfinite(external_load_direction_deg) || !std::isfinite(cfg.external_load_z)) {
        MessageBoxW(g_hwnd, L"external load values must be finite.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    if (external_load_magnitude < 0.0) { MessageBoxW(g_hwnd, L"external load magnitude cannot be negative.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    if (!std::isfinite(cfg.bearing_mass) ||
        (cfg.motion_model == MotionModel::MOVING_BEARING && cfg.bearing_mass <= 0.0)) {
        MessageBoxW(g_hwnd, L"bearing_mass must be positive when MOVING_BEARING is selected.", L"Pancake", MB_ICONERROR | MB_OK);
        return false;
    }
    for (double value : {cfg.bearing_initial_x, cfg.bearing_initial_y, cfg.bearing_initial_z,
                         cfg.bearing_initial_vx, cfg.bearing_initial_vy, cfg.bearing_initial_vz,
                         cfg.bearing_stiffness_x, cfg.bearing_stiffness_y, cfg.bearing_stiffness_z,
                         cfg.bearing_damping_x, cfg.bearing_damping_y, cfg.bearing_damping_z}) {
        if (!std::isfinite(value)) {
            MessageBoxW(g_hwnd, L"bearing motion support and initial-state values must be finite.", L"Pancake", MB_ICONERROR | MB_OK);
            return false;
        }
    }
    if (!std::isfinite(cfg.min_film_thickness) || cfg.min_film_thickness <= 0.0) {
        MessageBoxW(g_hwnd, L"min_film_thickness must be positive.", L"Pancake", MB_ICONERROR | MB_OK);
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
        << "write_interval = " << native_or_converted_text(kEditWriteInterval, kComboUnitWriteInterval, cfg.write_interval, unit_is_native_time(kComboUnitWriteInterval)) << "\n"
        << "# solution_mode options: TRANSIENT, STEADY_STATE\n"
        << "solution_mode = " << to_config_value(cfg.solution_mode) << "\n"
        << "# Options: EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICOLSON, RK2, RK4\n"
        << "pressure_time_method = " << to_config_value(cfg.pressure_time_method) << "\n"
        << "motion_time_method = " << to_config_value(cfg.motion_time_method) << "\n"
        << "temperature_time_method = " << to_config_value(cfg.temperature_time_method) << "\n\n";

    out << "# Physics\n"
        << "omega = " << (combo_is(kComboUnitOmega, L"rad/s") ? trim_copy(get_edit_utf8(kEditOmega)) : compact_number(cfg.omega)) << "\n"
        << "mu = " << trim_copy(get_edit_utf8(kEditMu)) << "\n"
        << "rho = " << trim_copy(get_edit_utf8(kEditRho)) << "\n"
        << "p_cav = " << native_or_converted_text(kEditPCav, kComboUnitPCav, cfg.p_cav, unit_is_native_pressure(kComboUnitPCav)) << "\n"
        << "bulk_modulus = " << native_or_converted_text(kEditBulkModulus, kComboUnitBulkModulus, cfg.bulk_modulus, unit_is_native_pressure(kComboUnitBulkModulus)) << "\n\n";

    out << "# Energy Equation\n"
        << "# temperature_model options: ISOTHERMAL, ENERGY_EQUATION\n"
        << "temperature_model = " << to_config_value(cfg.temperature_model) << "\n"
        << "temperature_initial = " << trim_copy(get_edit_utf8(kEditTemperatureInitial)) << "\n"
        << "temperature_reference = " << trim_copy(get_edit_utf8(kEditTemperatureReference)) << "\n"
        << "journal_wall_temperature = " << trim_copy(get_edit_utf8(kEditJournalWallTemperature)) << "\n"
        << "bearing_wall_temperature = " << trim_copy(get_edit_utf8(kEditBearingWallTemperature)) << "\n"
        << "rho_cp = " << trim_copy(get_edit_utf8(kEditRhoCp)) << "\n"
        << "thermal_conductivity = " << trim_copy(get_edit_utf8(kEditThermalConductivity)) << "\n"
        << "journal_heat_transfer = " << trim_copy(get_edit_utf8(kEditJournalHeatTransfer)) << "\n"
        << "bearing_heat_transfer = " << trim_copy(get_edit_utf8(kEditBearingHeatTransfer)) << "\n\n";

    out << "# Bearing Motion\n"
        << "# motion_model options: STATIC, MOVING_BEARING\n"
        << "motion_model = " << to_config_value(cfg.motion_model) << "\n"
        << "bearing_initial_from_attitude = " << (cfg.bearing_initial_from_attitude ? "true" : "false") << "\n"
        << "bearing_initial_x = " << compact_number(cfg.bearing_initial_x) << "\n"
        << "bearing_initial_y = " << compact_number(cfg.bearing_initial_y) << "\n"
        << "bearing_initial_z = " << compact_number(cfg.bearing_initial_z) << "\n"
        << "bearing_initial_vx = " << compact_number(cfg.bearing_initial_vx) << "\n"
        << "bearing_initial_vy = " << compact_number(cfg.bearing_initial_vy) << "\n"
        << "bearing_initial_vz = " << compact_number(cfg.bearing_initial_vz) << "\n"
        << "bearing_mass = " << compact_number(cfg.bearing_mass) << "\n"
        << "bearing_stiffness_x = " << compact_number(cfg.bearing_stiffness_x) << "\n"
        << "bearing_stiffness_y = " << compact_number(cfg.bearing_stiffness_y) << "\n"
        << "bearing_stiffness_z = " << compact_number(cfg.bearing_stiffness_z) << "\n"
        << "bearing_damping_x = " << compact_number(cfg.bearing_damping_x) << "\n"
        << "bearing_damping_y = " << compact_number(cfg.bearing_damping_y) << "\n"
        << "bearing_damping_z = " << compact_number(cfg.bearing_damping_z) << "\n"
        << "external_load_x = " << compact_number(cfg.external_load_x) << "\n"
        << "external_load_y = " << compact_number(cfg.external_load_y) << "\n"
        << "external_load_z = " << compact_number(cfg.external_load_z) << "\n"
        << "min_film_thickness = " << compact_number(cfg.min_film_thickness) << "\n"
        << "stop_on_nonpositive_film = " << (cfg.stop_on_nonpositive_film ? "true" : "false") << "\n\n";

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
    const std::string inlet_theta = combo_is(kComboUnitInletTheta, L"deg") ? trim_copy(get_edit_utf8(kEditInletTheta)) : compact_number(inlet.theta);
    const std::string inlet_z = unit_is_native_length(kComboUnitInletZ) ? trim_copy(get_edit_utf8(kEditInletZ)) : compact_number(inlet.z);
    const std::string inlet_size = (inlet.type == InletConfig::Type::CIRCULAR)
        ? native_or_converted_text(kEditInletSize, kComboUnitInletSize, inlet.size, unit_is_native_length(kComboUnitInletSize))
        : (combo_is(kComboUnitInletSize, L"deg") ? trim_copy(get_edit_utf8(kEditInletSize)) : compact_number(inlet.size));
    const std::string inlet_pressure = native_or_converted_text(kEditInletPressure, kComboUnitInletPressure, inlet.p_supply, unit_is_native_pressure(kComboUnitInletPressure));
    if (!cfg.inlets.empty() && cfg.inlets.front().type == InletConfig::Type::CIRCULAR) {
        out << "inlet_circular = " << inlet_theta << " " << inlet_z << " " << inlet_size << " " << inlet_pressure << "\n"
            << "# inlet_groove = " << inlet_theta << " " << compact_number(inlet.size) << " " << inlet_pressure << "\n";
    } else {
        out << "inlet_groove = " << inlet_theta << " " << inlet_size << " " << inlet_pressure << "\n"
            << "# inlet_circular = " << inlet_theta << " " << compact_number(0.5 * cfg.L) << " " << compact_number(cfg.c) << " " << inlet_pressure << "\n";
    }
    out << "\n";

    out << "# Output\n"
        << "# output_fields accepts comma-separated field names.\n"
        << "# Options: pressure, film_content, h, rho, temperature, heat_generation,\n"
        << "#          inlet_indicator, velocity,\n"
        << "#          pressure_force_x, pressure_force_y, pressure_force_z,\n"
        << "#          load_x, load_y, load_z (legacy pressure-force aliases),\n"
        << "#          viscous_force_x, viscous_force_y, viscous_force_z,\n"
        << "#          fluid_force_x, fluid_force_y, fluid_force_z,\n"
        << "#          external_load_x, external_load_y, external_load_z,\n"
        << "#          bearing_x, bearing_y, bearing_z, bearing_attitude_angle,\n"
        << "#          friction_torque\n"
        << "output_dir = " << trim_copy(get_edit_utf8(kEditOutputDir)) << "\n"
        << "filename_prefix = " << trim_copy(get_edit_utf8(kEditFilenamePrefix)) << "\n"
        << "output_write_3d = " << (cfg.output_write_3d ? "true" : "false") << "\n"
        << "output_write_flat = " << (cfg.output_write_flat ? "true" : "false") << "\n"
        << "output_fields = " << cfg.output_fields_text() << "\n";
    return out.str();
}

// ----- presets -------------------------------------------------------------
SimulationConfig preset_config(int index) {
    // Base on a clean aligned bearing (struct defaults) and vary one knob each.
    SimulationConfig cfg;            // R=0.01 c=1e-4 e=8e-5 L=0.02 ...
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0; groove.size = 10.0; groove.p_supply = 2.0e5;
    cfg.inlets = {groove};
    switch (index) {
        case 1: // High eccentricity
            cfg.e = 0.95 * cfg.c;
            break;
        case 2: // Misaligned (small axial tilt)
            cfg.tilt_slope_x = 5.0e-4;
            break;
        case 3: { // Circular oil-supply hole
            InletConfig hole;
            hole.type = InletConfig::Type::CIRCULAR;
            hole.theta = 90.0; hole.z = 0.5 * cfg.L; hole.size = cfg.c; hole.p_supply = 2.0e5;
            cfg.inlets = {hole};
            break;
        }
        default: break; // Default bearing
    }
    return cfg;
}

void load_preset(int index) {
    g_config = preset_config(index);
    g_config_value_text.clear();
    g_config_inlets_text.clear();
    populate_form_from_config(g_config);
    g_dirty = true;
    recompute_summary();
}

// ----- unit system switching ----------------------------------------------
void convert_length_unit(int edit_id, int unit_id, const wchar_t* new_unit) {
    double v = 0.0;
    if (try_double(edit_id, v)) {
        const double si = v * unit_scale_to_meters(unit_id);
        set_combo_text(unit_id, new_unit);
        set_edit_text(edit_id, compact_number(si / unit_scale_to_meters(unit_id)));
    } else {
        set_combo_text(unit_id, new_unit);
    }
}
void convert_pressure_unit(int edit_id, int unit_id, const wchar_t* new_unit) {
    double v = 0.0;
    if (try_double(edit_id, v)) {
        const double si = v * unit_scale_to_pascals(unit_id);
        set_combo_text(unit_id, new_unit);
        set_edit_text(edit_id, compact_number(si / unit_scale_to_pascals(unit_id)));
    } else {
        set_combo_text(unit_id, new_unit);
    }
}
void convert_omega_unit(const wchar_t* new_unit) {
    double v = 0.0;
    if (try_double(kEditOmega, v)) {
        const double si = combo_is(kComboUnitOmega, L"rpm") ? v * 2.0 * kPi / 60.0 : v;
        set_combo_text(kComboUnitOmega, new_unit);
        const bool rpm = (wcscmp(new_unit, L"rpm") == 0);
        set_edit_text(kEditOmega, compact_number(rpm ? si * 60.0 / (2.0 * kPi) : si));
    } else {
        set_combo_text(kComboUnitOmega, new_unit);
    }
}

// Applies the chosen global system to per-field override combos in place,
// preserving the physical value (only the displayed unit/number change).
void apply_unit_system(bool engineering) {
    g_loading = true;
    const wchar_t* len = engineering ? L"mm" : L"m";
    const wchar_t* prs = engineering ? L"MPa" : L"Pa";
    for (int id : {kEditR, kEditC, kEditE, kEditL}) {
        const int unit = kComboUnitR + (id - kEditR);
        convert_length_unit(id, unit, len);
    }
    convert_length_unit(kEditInletZ, kComboUnitInletZ, len);
    if (combo_is(kComboInletType, L"CIRCULAR")) convert_length_unit(kEditInletSize, kComboUnitInletSize, len);
    convert_pressure_unit(kEditPCav, kComboUnitPCav, prs);
    convert_pressure_unit(kEditBulkModulus, kComboUnitBulkModulus, prs);
    convert_pressure_unit(kEditBcSouthVal, kComboUnitBcSouthVal, prs);
    convert_pressure_unit(kEditBcNorthVal, kComboUnitBcNorthVal, prs);
    convert_pressure_unit(kEditInletPressure, kComboUnitInletPressure, prs);
    convert_omega_unit(engineering ? L"rpm" : L"rad/s");
    g_loading = false;
    recompute_summary();
}

// ----- live derived quantities + validation --------------------------------
std::wstring fmt_eng(double value, int sig = 4) {
    if (!std::isfinite(value)) return L"--";
    std::ostringstream out;
    const double m = std::abs(value);
    if (m != 0.0 && (m < 1.0e-3 || m >= 1.0e5)) out << std::setprecision(sig) << std::scientific << value;
    else out << std::setprecision(sig) << std::defaultfloat << value;
    return utf8_to_wide(out.str());
}

void recompute_summary() {
    g_summary_pairs.clear();
    g_problems.clear();
    g_invalid_ids.clear();
    g_has_errors = false;

    auto err = [&](int id, const std::wstring& msg) {
        if (id) g_invalid_ids.insert(id);
        g_problems.push_back(L"✗ " + msg);
        g_has_errors = true;
    };
    auto warn = [&](const std::wstring& msg) { g_problems.push_back(L"⚠ " + msg); };

    const double R = si_length(kEditR, kComboUnitR);
    const double c = si_length(kEditC, kComboUnitC);
    const double e = si_length(kEditE, kComboUnitE);
    const double L = si_length(kEditL, kComboUnitL);
    const double omega = si_omega();

    if (!std::isfinite(R) || R <= 0.0) err(kEditR, L"R must be a positive length");
    if (!std::isfinite(c) || c <= 0.0) err(kEditC, L"c (clearance) must be positive");
    if (!std::isfinite(e) || e < 0.0) err(kEditE, L"e must be ≥ 0");
    if (!std::isfinite(L) || L <= 0.0) err(kEditL, L"L must be positive");
    if (std::isfinite(c) && std::isfinite(e) && c > 0.0 && e >= c)
        err(kEditE, L"e ≥ c: eccentricity ratio ε must be < 1");

    int n_theta = 0, n_z = 0;
    if (!try_int(kEditNTheta, n_theta) || n_theta <= 0) err(kEditNTheta, L"n_theta must be a positive integer");
    if (!try_int(kEditNZ, n_z) || n_z <= 0) err(kEditNZ, L"n_z must be a positive integer");

    double end_t = 0, dt = 0, wint = 0;
    const bool ok_end = try_double(kEditEndT, end_t) && end_t > 0.0;
    const bool ok_dt = try_double(kEditDt, dt) && dt > 0.0;
    const bool ok_wi = try_double(kEditWriteInterval, wint) && wint > 0.0;
    if (!ok_end) err(kEditEndT, L"end_t must be positive");
    if (!ok_dt) err(kEditDt, L"dt must be positive");
    if (!ok_wi) err(kEditWriteInterval, L"write_interval must be positive");
    if (ok_dt && ok_end && dt > end_t) warn(L"dt > end_t: only one step will run");
    if (ok_dt && ok_wi && dt > wint) warn(L"dt > write_interval: outputs every step regardless");
    if (ok_wi && ok_end && wint > end_t) warn(L"write_interval > end_t: only the initial frame is written");

    double mu = 0, rho = 0;
    if (!try_double(kEditMu, mu) || mu <= 0.0) err(kEditMu, L"mu (viscosity) must be positive");
    if (!try_double(kEditRho, rho) || rho <= 0.0) err(kEditRho, L"rho (density) must be positive");

    const bool steady_state = combo_is(kComboSolutionMode, L"STEADY_STATE");
    const bool energy_equation = combo_is(kComboTemperatureModel, L"ENERGY_EQUATION");
    for (int id : {kEditTemperatureInitial, kEditTemperatureReference,
                   kEditJournalWallTemperature, kEditBearingWallTemperature}) {
        double value = 0.0;
        if (!try_double(id, value) || !std::isfinite(value)) err(id, L"temperature values must be numeric");
    }
    double rho_cp = 0.0, thermal_conductivity = 0.0;
    double journal_heat_transfer = 0.0, bearing_heat_transfer = 0.0;
    if (!try_double(kEditRhoCp, rho_cp) || !std::isfinite(rho_cp) || (energy_equation && rho_cp <= 0.0))
        err(kEditRhoCp, L"rho_cp must be positive for energy-equation runs");
    if (!try_double(kEditThermalConductivity, thermal_conductivity) ||
        !std::isfinite(thermal_conductivity) || (energy_equation && thermal_conductivity <= 0.0))
        err(kEditThermalConductivity, L"thermal_conductivity must be positive for energy-equation runs");
    if (!try_double(kEditJournalHeatTransfer, journal_heat_transfer) ||
        !std::isfinite(journal_heat_transfer) || journal_heat_transfer < 0.0)
        err(kEditJournalHeatTransfer, L"journal heat-transfer coefficient must be >= 0");
    if (!try_double(kEditBearingHeatTransfer, bearing_heat_transfer) ||
        !std::isfinite(bearing_heat_transfer) || bearing_heat_transfer < 0.0)
        err(kEditBearingHeatTransfer, L"bearing heat-transfer coefficient must be >= 0");
    if (steady_state && energy_equation && journal_heat_transfer + bearing_heat_transfer <= 0.0)
        warn(L"steady energy solve has no wall heat sink; add heat-transfer coefficients");

    double external_load_magnitude = 0.0, external_load_direction = 0.0, external_load_z = 0.0;
    const bool ok_load_mag = try_double(kEditExternalLoadMagnitude, external_load_magnitude) &&
        std::isfinite(external_load_magnitude) && external_load_magnitude >= 0.0;
    const bool ok_load_dir = try_double(kEditExternalLoadDirection, external_load_direction) && std::isfinite(external_load_direction);
    const bool ok_load_z = try_double(kEditExternalLoadZ, external_load_z) && std::isfinite(external_load_z);
    if (!ok_load_mag) err(kEditExternalLoadMagnitude, L"external load magnitude must be >= 0");
    if (!ok_load_dir) err(kEditExternalLoadDirection, L"external load direction must be numeric");
    if (!ok_load_z) err(kEditExternalLoadZ, L"external load z must be numeric");

    const bool moving_bearing = combo_is(kComboMotionModel, L"MOVING_BEARING");
    double bearing_mass = 0.0, min_film_thickness = 0.0;
    if (!try_double(kEditBearingMass, bearing_mass) || !std::isfinite(bearing_mass) || (moving_bearing && bearing_mass <= 0.0))
        err(kEditBearingMass, L"bearing mass must be positive for moving-bearing runs");
    if (!try_double(kEditMinFilmThickness, min_film_thickness) || !std::isfinite(min_film_thickness) || min_film_thickness <= 0.0)
        err(kEditMinFilmThickness, L"minimum film thickness must be positive");
    for (int id : {kEditBearingInitialX, kEditBearingInitialY, kEditBearingInitialZ,
                   kEditBearingInitialVx, kEditBearingInitialVy, kEditBearingInitialVz,
                   kEditBearingStiffnessX, kEditBearingStiffnessY, kEditBearingStiffnessZ,
                   kEditBearingDampingX, kEditBearingDampingY, kEditBearingDampingZ}) {
        double value = 0.0;
        if (!try_double(id, value) || !std::isfinite(value)) err(id, L"motion support and initial-state values must be numeric");
    }

    if (trim_copy(get_edit_utf8(kEditOutputDir)).empty()) err(kEditOutputDir, L"output_dir cannot be empty");
    if (trim_copy(get_edit_utf8(kEditFilenamePrefix)).empty()) err(kEditFilenamePrefix, L"filename_prefix cannot be empty");

    // Derived quantities (only when the underlying inputs are sane).
    if (std::isfinite(c) && c > 0.0 && std::isfinite(e) && e >= 0.0) {
        const double eps = e / c;
        g_summary_pairs.push_back({L"ε = e/c", fmt_eng(eps) + (eps < 1.0 ? L"" : L"  (!)")});
        g_summary_pairs.push_back({L"h_min", fmt_eng((c - e) * 1.0e6) + L" µm"});
        g_summary_pairs.push_back({L"h_max", fmt_eng((c + e) * 1.0e6) + L" µm"});
    }
    if (std::isfinite(c) && std::isfinite(R) && R > 0.0)
        g_summary_pairs.push_back({L"c/R", fmt_eng(c / R)});
    if (std::isfinite(L) && std::isfinite(R) && R > 0.0)
        g_summary_pairs.push_back({L"L/D", fmt_eng(L / (2.0 * R))});
    if (std::isfinite(omega) && std::isfinite(R)) {
        g_summary_pairs.push_back({L"U = ωR", fmt_eng(omega * R) + L" m/s"});
        g_summary_pairs.push_back({L"speed", fmt_eng(omega * 60.0 / (2.0 * kPi)) + L" rpm"});
    }
    g_summary_pairs.push_back({L"solution", steady_state ? L"steady state" : L"transient"});
    g_summary_pairs.push_back({L"thermal", energy_equation ? L"energy equation" : L"isothermal"});
    if (ok_load_mag && ok_load_dir && ok_load_z) {
        const double load_rad = external_load_direction * kPi / 180.0;
        g_summary_pairs.push_back({L"|F_ext,xy|", fmt_eng(external_load_magnitude) + L" N"});
        g_summary_pairs.push_back({L"F_ext angle", fmt_eng(external_load_direction) + L" deg"});
        g_summary_pairs.push_back({L"F_ext z", fmt_eng(external_load_z) + L" N"});
        g_summary_pairs.push_back({L"F_ext x/y", fmt_eng(external_load_magnitude * std::cos(load_rad)) + L" / " +
                                      fmt_eng(external_load_magnitude * std::sin(load_rad)) + L" N"});
    }
    if (g_preview.ok) {
        g_summary_pairs.push_back({utf8_to_wide(g_preview.field) + L" mean", fmt_eng(g_preview.mean_value)});
        g_summary_pairs.push_back({utf8_to_wide(g_preview.field) + L" max", fmt_eng(g_preview.max_value)});
    }

    if (g_summary != nullptr) InvalidateRect(g_summary, nullptr, TRUE);
    if (g_run_button != nullptr) {
        const bool can_run = !g_has_errors && !g_solver.running;
        EnableWindow(g_run_button, can_run ? TRUE : FALSE);
    }
    // Repaint invalid-marked edits.
    for (const auto& ref : g_control_refs) {
        if (ref.id >= kEditR && ref.id <= kEditLastParameter) InvalidateRect(ref.hwnd, nullptr, TRUE);
    }
}

// ----- dirty / title / save / load ----------------------------------------
void update_title() {
    std::wstring title = L"Pancake";
    if (!g_config_path.empty()) title += L" — " + g_config_path.filename().wstring();
    if (g_dirty) title += L" *";
    if (g_hwnd != nullptr) SetWindowTextW(g_hwnd, title.c_str());
}

void mark_dirty() {
    if (g_loading) return;
    g_dirty = true;
    update_title();
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

bool g_active_raw = false;   // whether the Raw editor currently drives Save

bool save_to_path(const fs::path& path) {
    SimulationConfig cfg = g_config;
    std::string config_text;
    if (g_active_raw) {
        config_text = get_edit_utf8(kEditRawConfig);
        cfg = SimulationConfig{};
        try { cfg.load_from_text(config_text); }
        catch (const std::exception& ex) {
            MessageBoxW(g_hwnd, utf8_to_wide(ex.what()).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
            return false;
        }
    } else if (!apply_form_to_config(cfg)) {
        return false;
    } else {
        config_text = form_config_text(cfg);
    }
    if (!write_text_file(path, config_text)) return false;

    g_config = cfg;
    g_config_path = path;
    g_config_value_text = parse_config_value_text(config_text);
    g_config_inlets_text = parse_config_inlets_text(config_text);
    populate_form_from_config(g_config);
    recompute_summary();
    g_dirty = false;
    update_title();
    return true;
}

bool save_config() {
    if (g_config_path.empty()) g_config_path = g_executable_dir / L"config.txt";
    return save_to_path(g_config_path);
}

std::wstring run_file_dialog(bool save, const wchar_t* filter, const wchar_t* default_ext, const fs::path& initial) {
    wchar_t file[MAX_PATH] = {};
    if (!initial.empty()) {
        const std::wstring s = initial.wstring();
        wcsncpy(file, s.c_str(), MAX_PATH - 1);
    }
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_hwnd;
    ofn.lpstrFilter = filter;
    ofn.lpstrFile = file;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = default_ext;
    ofn.Flags = OFN_NOCHANGEDIR | OFN_PATHMUSTEXIST | (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
    const BOOL ok = save ? GetSaveFileNameW(&ofn) : GetOpenFileNameW(&ofn);
    return ok ? std::wstring(file) : std::wstring();
}

bool confirm_discard_changes() {
    if (!g_dirty) return true;
    const int result = MessageBoxW(g_hwnd, L"Discard the unsaved configuration changes?", L"Pancake",
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
    }
    populate_form_from_config(g_config);
    recompute_summary();
    g_dirty = false;
    update_title();
}

void open_config_dialog() {
    if (!confirm_discard_changes()) return;
    const std::wstring chosen = run_file_dialog(false,
        L"Config files (*.txt)\0*.txt\0All files (*.*)\0*.*\0", L"txt", g_config_path);
    if (chosen.empty()) return;
    g_config_path = chosen;
    load_config();
}

void save_as_dialog() {
    const std::wstring chosen = run_file_dialog(true,
        L"Config files (*.txt)\0*.txt\0All files (*.*)\0*.*\0", L"txt",
        g_config_path.empty() ? (g_executable_dir / L"config.txt") : g_config_path);
    if (chosen.empty()) return;
    save_to_path(chosen);
}

void reset_defaults() {
    if (!confirm_discard_changes()) return;
    g_config = SimulationConfig{};
    g_config_value_text.clear();
    g_config_inlets_text.clear();
    populate_form_from_config(g_config);
    recompute_summary();
    g_dirty = true;
    update_title();
}

void browse_solver() {
    const std::wstring chosen = run_file_dialog(false,
        L"Executable (*.exe)\0*.exe\0All files (*.*)\0*.*\0", L"exe",
        g_solver_path.empty() ? (g_executable_dir / L"pancake.exe") : fs::path(g_solver_path));
    if (!chosen.empty()) g_solver_path = chosen;
}

fs::path output_directory_for_config(const SimulationConfig& cfg) {
    fs::path output_path = utf8_to_wide(cfg.output_dir);
    if (output_path.is_relative()) output_path = g_executable_dir / output_path;
    return output_path;
}

void refresh_preview(bool reload_steps = false);

void sync_raw_from_form() {
    SimulationConfig cfg = g_config;
    if (!apply_form_to_config(cfg)) return;
    g_config = cfg;
    g_loading = true;
    set_edit_text(kEditRawConfig, form_config_text(g_config));
    g_loading = false;
}

void apply_raw_to_form() {
    const std::string raw_text = get_edit_utf8(kEditRawConfig);
    SimulationConfig cfg;
    try { cfg.load_from_text(raw_text); }
    catch (const std::exception& ex) {
        MessageBoxW(g_hwnd, utf8_to_wide(ex.what()).c_str(), L"Pancake", MB_ICONERROR | MB_OK);
        return;
    }
    g_config = cfg;
    g_config_value_text = parse_config_value_text(raw_text);
    g_config_inlets_text = parse_config_inlets_text(raw_text);
    populate_form_from_config(g_config);
    recompute_summary();
    mark_dirty();
}

// ----- log (ANSI-coloured) -------------------------------------------------
bool try_parse_ansi_codes(const std::wstring& text, size_t offset, size_t& next, std::vector<int>& codes) {
    size_t start = std::wstring::npos;
    if (text[offset] == L'\x1b' && offset + 1 < text.size() && text[offset + 1] == L'[') start = offset + 2;
    else if (text[offset] == L'[' && offset + 2 < text.size() && std::iswdigit(text[offset + 1])) start = offset + 1;
    else return false;

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
        if (!token.empty()) { try { codes.push_back(std::stoi(token)); } catch (...) { return false; } }
        if (code_end == std::wstring::npos) break;
        code_start = code_end + 1;
    }
    next = end + 1;
    return !codes.empty();
}

void apply_ansi_codes(const std::vector<int>& codes, COLORREF& color) {
    for (int code : codes) {
        switch (code) {
            case 0: case 39: color = RGB(38, 46, 55); break;
            case 31: case 91: color = RGB(183, 61, 69); break;
            case 32: case 92: color = RGB(25, 128, 74); break;
            case 33: case 93: color = RGB(153, 109, 21); break;
            case 34: case 94: color = RGB(40, 110, 196); break;
            case 35: case 95: color = RGB(150, 70, 168); break;
            case 36: case 96: color = RGB(20, 119, 145); break;
            case 90: color = RGB(99, 110, 123); break;
            default: break;
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

// Pull the latest "t=<value>" progress token out of solver output.
void update_progress_from_text(const std::wstring& text) {
    size_t pos = 0;
    double latest = -1.0;
    while ((pos = text.find(L"t=", pos)) != std::wstring::npos) {
        size_t j = pos + 2;
        std::wstring num;
        while (j < text.size() && (std::iswdigit(text[j]) || text[j] == L'.' || text[j] == L'e' ||
               text[j] == L'E' || text[j] == L'+' || text[j] == L'-')) {
            num += text[j++];
        }
        try { latest = std::stod(num); } catch (...) {}
        pos = j;
    }
    if (latest >= 0.0 && g_config.end_t > 0.0) {
        g_run_progress = std::clamp(latest / g_config.end_t, 0.0, 1.0);
    }
    if (text.find(L"finished successfully") != std::wstring::npos) g_run_progress = 1.0;
}

// ----- solver process ------------------------------------------------------
std::wstring quote_arg(const std::wstring& arg) {
    std::wstring result = L"\"";
    for (wchar_t ch : arg) { if (ch == L'"') result += L'\\'; result += ch; }
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
    for (const auto& candidate : candidates) if (fs::exists(candidate)) return candidate;
    wchar_t path[MAX_PATH] = {};
    const DWORD length = SearchPathW(nullptr, L"mpiexec.exe", nullptr, MAX_PATH, path, nullptr);
    if (length > 0 && length < MAX_PATH) return fs::path(path);
    return {};
}

void close_solver_handles() {
    if (g_solver.stdout_read != nullptr) { CloseHandle(g_solver.stdout_read); g_solver.stdout_read = nullptr; }
    if (g_solver.process_info.hThread != nullptr) { CloseHandle(g_solver.process_info.hThread); g_solver.process_info.hThread = nullptr; }
    if (g_solver.process_info.hProcess != nullptr) { CloseHandle(g_solver.process_info.hProcess); g_solver.process_info.hProcess = nullptr; }
    g_solver.process_info.dwProcessId = 0;
    g_solver.running = false;
}

void invalidate_action_bar() {
    if (g_hwnd == nullptr) return;
    RECT r{};
    GetClientRect(g_hwnd, &r);
    r.bottom = kActionBarH;
    InvalidateRect(g_hwnd, &r, TRUE);
}

void set_run_state(RunState state) {
    g_run_state = state;
    const bool running = (state == RunState::Running);
    EnableWindow(g_run_button, (!running && !g_has_errors) ? TRUE : FALSE);
    EnableWindow(g_stop_button, running ? TRUE : FALSE);
    if (g_progress != nullptr) {
        ShowWindow(g_progress, running ? SW_SHOW : SW_HIDE);
        if (running) SendMessageW(g_progress, PBM_SETPOS, static_cast<WPARAM>(g_run_progress * 100.0), 0);
    }
    invalidate_action_bar();
}

bool read_mpi_ranks(int& ranks) {
    if (!read_int(kEditMpiRanks, L"MPI ranks", ranks)) return false;
    if (ranks < 1) { MessageBoxW(g_hwnd, L"MPI ranks must be at least 1.", L"Pancake", MB_ICONERROR | MB_OK); return false; }
    return true;
}

fs::path project_root_from_gui() {
    const fs::path parent = g_executable_dir.parent_path();
    if (fs::exists(parent / L"CMakeLists.txt")) return parent;
    return fs::current_path();
}

bool make_solver_command(int ranks, std::wstring& command, fs::path& working_dir, std::wstring& missing_message) {
    const fs::path project_root = project_root_from_gui();
    std::vector<fs::path> native_candidates;
    if (!g_solver_path.empty()) {
        fs::path requested = g_solver_path;
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
        L"The GUI expects pancake.exe beside pancake_gui.exe or under build-windows-mingw.\n"
        L"Use \"Solver...\" to point at a custom build.";
    return false;
}

void start_solver() {
    if (g_solver.running) return;
    if (!save_config()) return;
    if (g_has_errors) { MessageBoxW(g_hwnd, L"Fix the highlighted parameters before running.", L"Pancake", MB_ICONWARNING | MB_OK); return; }

    int ranks = 1;
    if (!read_mpi_ranks(ranks)) return;

    std::wstring command, missing_message;
    fs::path working_dir = g_executable_dir;
    if (!make_solver_command(ranks, command, working_dir, missing_message)) {
        MessageBoxW(g_hwnd, missing_message.c_str(), L"Pancake", MB_ICONINFORMATION | MB_OK);
        return;
    }

    SECURITY_ATTRIBUTES security{};
    security.nLength = sizeof(security);
    security.bInheritHandle = TRUE;
    HANDLE read_pipe = nullptr, write_pipe = nullptr;
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
    const BOOL started = CreateProcessW(nullptr, mutable_command.data(), nullptr, nullptr, TRUE,
                                        CREATE_NO_WINDOW, nullptr, working_dir.c_str(), &startup, &process_info);
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
    g_run_progress = 0.0;
    g_run_elapsed_s = 0.0;
    g_run_start_tick = GetTickCount64();
    set_run_state(RunState::Running);
    SetTimer(g_hwnd, kProcessTimer, 400, nullptr);
}

void stop_solver() {
    if (!g_solver.running) return;
    TerminateProcess(g_solver.process_info.hProcess, 1);
    append_log(L"\r\nSolver stop requested.\r\n");
    close_solver_handles();
    KillTimer(g_hwnd, kProcessTimer);
    set_run_state(RunState::Idle);
}

void poll_solver() {
    if (!g_solver.running) return;
    g_run_elapsed_s = (GetTickCount64() - g_run_start_tick) / 1000.0;

    char buffer[4096];
    DWORD available = 0;
    while (g_solver.stdout_read != nullptr &&
           PeekNamedPipe(g_solver.stdout_read, nullptr, 0, nullptr, &available, nullptr) && available > 0) {
        DWORD bytes_read = 0;
        const DWORD to_read = std::min<DWORD>(available, sizeof(buffer) - 1);
        if (!ReadFile(g_solver.stdout_read, buffer, to_read, &bytes_read, nullptr) || bytes_read == 0) break;
        buffer[bytes_read] = '\0';
        const std::wstring chunk = utf8_to_wide(std::string(buffer, bytes_read));
        append_log(chunk);
        update_progress_from_text(chunk);
    }

    if (g_progress != nullptr) SendMessageW(g_progress, PBM_SETPOS, static_cast<WPARAM>(g_run_progress * 100.0), 0);
    invalidate_action_bar();

    DWORD exit_code = STILL_ACTIVE;
    if (GetExitCodeProcess(g_solver.process_info.hProcess, &exit_code) && exit_code != STILL_ACTIVE) {
        append_log(L"\r\nSolver exited with code " + std::to_wstring(exit_code) + L".\r\n");
        close_solver_handles();
        KillTimer(g_hwnd, kProcessTimer);
        g_run_progress = exit_code == 0 ? 1.0 : g_run_progress;
        set_run_state(exit_code == 0 ? RunState::Done : RunState::Failed);
        refresh_preview(true);
        recompute_summary();
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
    if (g_config.output_write_flat && fs::exists(output_dir / L"flat" / L"results.pvd")) pvd = output_dir / L"flat" / L"results.pvd";
    else pvd = output_dir / L"results.pvd";
    if (!fs::exists(pvd)) { MessageBoxW(g_hwnd, L"No PVD result file exists yet.", L"Pancake", MB_ICONINFORMATION | MB_OK); return; }
    ShellExecuteW(g_hwnd, L"open", pvd.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
}

// ----- result preview (VTS parsing + contour) ------------------------------
void populate_preview_fields() {
    HWND combo = find_control(kComboPreviewField);
    if (combo == nullptr) return;
    SendMessageW(combo, CB_RESETCONTENT, 0, 0);
    for (const auto& field : kOutputFields) {
        if (std::string(field.key) == "velocity") continue;
        SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(utf8_to_wide(field.key).c_str()));
    }
    SendMessageW(combo, CB_SETCURSEL, 0, 0);
}

bool parse_extent_values(const std::string& text, size_t extent_pos, int& x0, int& x1, int& y0, int& y1) {
    const size_t value_start = extent_pos + 8;
    const size_t value_end = text.find('"', value_start);
    if (value_end == std::string::npos) return false;
    std::stringstream values(text.substr(value_start, value_end - value_start));
    int z0 = 0, z1 = 0;
    if (!(values >> x0 >> x1 >> y0 >> y1 >> z0 >> z1)) return false;
    return x1 > x0 && y1 > y0;
}

bool parse_vts_extents(const std::string& text, int& whole_nx, int& whole_ny, int& piece_x0, int& piece_x1, int& piece_y0, int& piece_y1) {
    int wx0 = 0, wx1 = 0, wy0 = 0, wy1 = 0;
    const size_t whole_pos = text.find("WholeExtent=\"");
    if (whole_pos == std::string::npos || !parse_extent_values(text, whole_pos + 5, wx0, wx1, wy0, wy1)) return false;
    const size_t piece_tag = text.find("<Piece");
    if (piece_tag == std::string::npos) return false;
    const size_t piece_pos = text.find("Extent=\"", piece_tag);
    if (piece_pos == std::string::npos || !parse_extent_values(text, piece_pos, piece_x0, piece_x1, piece_y0, piece_y1)) return false;
    whole_nx = wx1 - wx0;
    whole_ny = wy1 - wy0;
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
    try { step = std::stoi(name.substr(step_start, step_end - step_start)); return step >= 0; }
    catch (...) { return false; }
}

bool is_processor_directory(const fs::path& path) {
    const std::wstring name = path.filename().wstring();
    return name.rfind(L"processor", 0) == 0;
}

void collect_preview_steps_from_root(const fs::path& root, const std::string& prefix, std::vector<PreviewStep>& steps) {
    if (!fs::exists(root)) return;
    std::map<int, PreviewStep> by_step;
    for (const auto& dir_entry : fs::directory_iterator(root)) {
        if (!dir_entry.is_directory() || !is_processor_directory(dir_entry.path())) continue;
        for (const auto& file_entry : fs::directory_iterator(dir_entry.path())) {
            if (!file_entry.is_regular_file()) continue;
            if (!file_name_matches_vts(file_entry.path(), prefix)) continue;
            int step = -1;
            if (!parse_step_from_vts_name(file_entry.path(), prefix, step)) continue;
            PreviewStep& item = by_step[step];
            item.step = step;
            item.root = root;
            item.files.push_back(file_entry.path());
            const auto write_time = file_entry.last_write_time();
            if (item.files.size() == 1 || write_time > item.write_time) item.write_time = write_time;
        }
    }
    for (auto& pair : by_step) {
        std::sort(pair.second.files.begin(), pair.second.files.end());
        steps.push_back(std::move(pair.second));
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
    collect_preview_steps_from_root(output_dir / L"flat", g_config.filename_prefix, g_preview_steps);
    if (g_preview_steps.empty()) collect_preview_steps_from_root(output_dir, g_config.filename_prefix, g_preview_steps);
    std::sort(g_preview_steps.begin(), g_preview_steps.end(), [](const PreviewStep& a, const PreviewStep& b) { return a.step < b.step; });

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

bool selected_preview_files(const PreviewStep*& selection, int& step) {
    if (g_preview_steps.empty()) return false;
    step = selected_preview_step();
    if (step < 0) {
        const auto latest = std::max_element(g_preview_steps.begin(), g_preview_steps.end(), [](const PreviewStep& a, const PreviewStep& b) { return a.step < b.step; });
        if (latest == g_preview_steps.end()) return false;
        step = latest->step;
        selection = &*latest;
        return true;
    }
    const auto item = std::find_if(g_preview_steps.begin(), g_preview_steps.end(), [step](const PreviewStep& c) { return c.step == step; });
    if (item == g_preview_steps.end()) return false;
    selection = &*item;
    return true;
}

bool load_preview_values_for_step(const PreviewStep& preview_step, int step, const std::string& field,
                                  int& nx, int& ny, std::vector<double>& values, std::wstring& source_label) {
    bool loaded_any = false;
    nx = 0; ny = 0; values.clear();
    for (const auto& path : preview_step.files) {
        const std::string text = read_file_text(path);
        int whole_nx = 0, whole_ny = 0, px0 = 0, px1 = 0, py0 = 0, py1 = 0;
        std::vector<double> local_values;
        if (!parse_vts_extents(text, whole_nx, whole_ny, px0, px1, py0, py1) || !parse_data_array(text, field, local_values)) continue;
        if (!loaded_any) {
            nx = whole_nx; ny = whole_ny;
            values.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), std::numeric_limits<double>::quiet_NaN());
            loaded_any = true;
        }
        if (whole_nx != nx || whole_ny != ny) continue;
        const int local_nx = px1 - px0;
        const int local_ny = py1 - py0;
        const size_t expected = static_cast<size_t>(local_nx) * static_cast<size_t>(local_ny);
        if (local_values.size() < expected) continue;
        for (int j = 0; j < local_ny; ++j) {
            const int gj = py0 + j;
            if (gj < 0 || gj >= ny) continue;
            for (int i = 0; i < local_nx; ++i) {
                const int gi = px0 + i;
                if (gi < 0 || gi >= nx) continue;
                values[static_cast<size_t>(gj) * nx + gi] = local_values[static_cast<size_t>(j) * local_nx + i];
            }
        }
    }
    if (!loaded_any) return false;
    source_label = preview_step.root.wstring() + L" (step " + std::to_wstring(step) +
                   L", " + std::to_wstring(preview_step.files.size()) + L" files)";
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
    refresh_preview(false);
}

COLORREF heat_color(double t) {
    t = std::clamp(t, 0.0, 1.0);
    struct Stop { double at; int r; int g; int b; };
    constexpr std::array<Stop, 4> stops = {{ {0.0, 33, 75, 99}, {0.36, 35, 151, 139}, {0.68, 241, 196, 83}, {1.0, 188, 64, 56} }};
    Stop a = stops.front(), b = stops.back();
    for (size_t i = 1; i < stops.size(); ++i) { if (t <= stops[i].at) { a = stops[i - 1]; b = stops[i]; break; } }
    const double local = (b.at > a.at) ? (t - a.at) / (b.at - a.at) : 0.0;
    return RGB(static_cast<int>(a.r + (b.r - a.r) * local), static_cast<int>(a.g + (b.g - a.g) * local), static_cast<int>(a.b + (b.b - a.b) * local));
}

std::string plot_number(double value) {
    if (!std::isfinite(value)) return "nan";
    const double magnitude = std::abs(value);
    std::ostringstream out;
    if (magnitude > 0.0 && (magnitude < 1.0e-3 || magnitude >= 1.0e5)) out << std::setprecision(3) << std::scientific << value;
    else out << std::setprecision(5) << std::defaultfloat << value;
    return out.str();
}
std::wstring plot_number_wide(double value) { return utf8_to_wide(plot_number(value)); }

std::wstring field_unit_label(const std::string& field) {
    if (field == "pressure") return L"Pa";
    if (field == "h") return L"m";
    if (field == "rho") return L"kg/m^3";
    if (field == "temperature") return L"K";
    if (field == "heat_generation") return L"W/m^2";
    if (field == "pressure_force_x" || field == "pressure_force_y" || field == "pressure_force_z" ||
        field == "load_x" || field == "load_y" || field == "load_z" ||
        field == "viscous_force_x" || field == "viscous_force_y" || field == "viscous_force_z" ||
        field == "fluid_force_x" || field == "fluid_force_y" || field == "fluid_force_z" ||
        field == "external_load_x" || field == "external_load_y" || field == "external_load_z") return L"N";
    if (field == "bearing_x" || field == "bearing_y" || field == "bearing_z") return L"m";
    if (field == "bearing_attitude_angle") return L"deg";
    if (field == "friction_torque") return L"N m";
    if (field == "velocity") return L"m/s";
    if (field == "film_content" || field == "inlet_indicator") return L"-";
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

void refresh_preview(bool reload_steps) {
    const std::wstring field_wide = get_combo_text(kComboPreviewField);
    const std::string field = normalise_key(wide_to_utf8(field_wide));
    if (field.empty()) return;

    if (reload_steps || g_preview_steps.empty()) refresh_preview_step_items();
    const PreviewStep* selected_files = nullptr;
    int selected_step = -1;
    if (!selected_preview_files(selected_files, selected_step) || selected_files == nullptr) {
        g_preview = PreviewData{};
        g_preview.message = L"No VTS results found yet. Run the solver.";
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }

    int nx = 0, ny = 0;
    std::vector<double> values;
    std::wstring source_label;
    if (!load_preview_values_for_step(*selected_files, selected_step, field, nx, ny, values, source_label)) {
        g_preview = PreviewData{};
        g_preview.message = L"Selected field was not found in the selected timestep.";
        g_preview.source_path = selected_files->root.wstring();
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }
    const size_t expected = static_cast<size_t>(nx) * static_cast<size_t>(ny);
    if (values.size() < expected) {
        g_preview = PreviewData{};
        g_preview.message = L"Latest VTS file has incomplete field data.";
        InvalidateRect(g_preview_canvas, nullptr, TRUE);
        return;
    }
    values.resize(expected);

    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    int finite_count = 0;
    for (double value : values) {
        if (!std::isfinite(value)) continue;
        min_value = std::min(min_value, value);
        max_value = std::max(max_value, value);
        sum += value;
        ++finite_count;
    }
    if (!std::isfinite(min_value) || !std::isfinite(max_value)) { min_value = 0.0; max_value = 0.0; }

    g_preview.ok = true;
    g_preview.nx = nx;
    g_preview.ny = ny;
    g_preview.step = selected_step;
    g_preview.min_value = min_value;
    g_preview.max_value = max_value;
    g_preview.mean_value = finite_count > 0 ? sum / finite_count : 0.0;
    g_preview.field = field;
    g_preview.values = std::move(values);
    g_preview.source_path = source_label;
    InvalidateRect(g_preview_canvas, nullptr, TRUE);
    recompute_summary();
}

void draw_preview(HWND hwnd, HDC hdc) {
    RECT rect{};
    GetClientRect(hwnd, &rect);
    FillRect(hdc, &rect, g_brush_card);
    SetBkMode(hdc, TRANSPARENT);
    SelectObject(hdc, g_font);
    SetTextColor(hdc, kColText);

    if (!g_preview.ok) {
        DrawTextW(hdc, g_preview.message.c_str(), -1, &rect, DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
        return;
    }

    RECT plot = rect;
    plot.left += 64; plot.top += 78; plot.right -= 178; plot.bottom -= 78;
    if (plot.right <= plot.left || plot.bottom <= plot.top) return;

    std::wstring title = utf8_to_wide(g_preview.field) + L"  \x2022  step " + std::to_wstring(g_preview.step);
    RECT title_rect{rect.left + 16, rect.top + 10, rect.right - 16, rect.top + 32};
    SelectObject(hdc, g_heading_font);
    DrawTextW(hdc, title.c_str(), -1, &title_rect, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    SelectObject(hdc, g_font);

    std::ostringstream stats;
    stats << "min " << plot_number(g_preview.min_value) << "      max " << plot_number(g_preview.max_value)
          << "      grid " << g_preview.nx << " x " << g_preview.ny;
    RECT stats_rect{rect.left + 16, rect.top + 40, rect.right - 16, rect.top + 60};
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
    MoveToEx(hdc, plot.left, plot.bottom, nullptr); LineTo(hdc, plot.right, plot.bottom);
    MoveToEx(hdc, plot.left, plot.top, nullptr); LineTo(hdc, plot.left, plot.bottom);
    SetTextColor(hdc, kColMuted);
    for (int tick = 0; tick <= 4; ++tick) {
        const int x = plot.left + tick * (plot.right - plot.left) / 4;
        MoveToEx(hdc, x, plot.bottom, nullptr); LineTo(hdc, x, plot.bottom + 5);
        RECT tick_rect{x - 24, plot.bottom + 6, x + 24, plot.bottom + 24};
        DrawTextW(hdc, std::to_wstring(tick * 90).c_str(), -1, &tick_rect, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
    }
    for (int tick = 0; tick <= 4; ++tick) {
        const int y = plot.bottom - tick * (plot.bottom - plot.top) / 4;
        MoveToEx(hdc, plot.left - 5, y, nullptr); LineTo(hdc, plot.left, y);
        const double z_value = g_config.L * static_cast<double>(tick) / 4.0;
        RECT tick_rect{rect.left + 4, y - 10, plot.left - 8, y + 10};
        DrawTextW(hdc, plot_number_wide(z_value).c_str(), -1, &tick_rect, DT_RIGHT | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
    }
    RECT x_label{plot.left, plot.bottom + 26, plot.right, plot.bottom + 44};
    DrawTextW(hdc, L"angle from +Y load reference (deg)", -1, &x_label, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
    HFONT y_font = CreateFontW(-16, 0, 900, 900, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    HGDIOBJ old_y_font = SelectObject(hdc, y_font);
    const UINT old_align = SetTextAlign(hdc, TA_CENTER | TA_BASELINE);
    TextOutW(hdc, rect.left + 18, (plot.top + plot.bottom) / 2, L"z [m]", 5);
    SetTextAlign(hdc, old_align);
    SelectObject(hdc, old_y_font);
    DeleteObject(y_font);
    SelectObject(hdc, old_pen);
    DeleteObject(axis_pen);

    RECT colorbar{};
    colorbar.left = plot.right + 24; colorbar.right = colorbar.left + 18; colorbar.top = plot.top; colorbar.bottom = plot.bottom;
    const int colorbar_height = std::max(1, static_cast<int>(colorbar.bottom - colorbar.top));
    for (int y = colorbar.top; y < colorbar.bottom; ++y) {
        const double t = 1.0 - static_cast<double>(y - colorbar.top) / colorbar_height;
        HPEN cp = CreatePen(PS_SOLID, 1, heat_color(t));
        HGDIOBJ old_cp = SelectObject(hdc, cp);
        MoveToEx(hdc, colorbar.left, y, nullptr); LineTo(hdc, colorbar.right, y);
        SelectObject(hdc, old_cp); DeleteObject(cp);
    }
    HPEN cb_border = CreatePen(PS_SOLID, 1, kColBorder);
    HGDIOBJ old_cb_pen = SelectObject(hdc, cb_border);
    HGDIOBJ old_cb_brush = SelectObject(hdc, GetStockObject(NULL_BRUSH));
    Rectangle(hdc, colorbar.left, colorbar.top, colorbar.right, colorbar.bottom);
    SelectObject(hdc, old_cb_brush); SelectObject(hdc, old_cb_pen); DeleteObject(cb_border);

    RECT color_title{plot.right + 8, colorbar.top - 24, rect.right - 8, colorbar.top - 4};
    SelectObject(hdc, g_heading_font);
    SetTextColor(hdc, kColText);
    DrawTextW(hdc, colorbar_title().c_str(), -1, &color_title, DT_CENTER | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    SelectObject(hdc, g_small_font);
    HPEN tick_pen = CreatePen(PS_SOLID, 1, RGB(120, 130, 143));
    HGDIOBJ old_tick_pen = SelectObject(hdc, tick_pen);
    const int n_color_ticks = 9;
    for (int tick = 0; tick < n_color_ticks; ++tick) {
        const double frac = (n_color_ticks == 1) ? 0.0 : static_cast<double>(tick) / (n_color_ticks - 1);
        const int y = colorbar.bottom - static_cast<int>(std::llround(frac * (colorbar.bottom - colorbar.top)));
        const double value = g_preview.min_value + frac * (g_preview.max_value - g_preview.min_value);
        MoveToEx(hdc, colorbar.right, y, nullptr); LineTo(hdc, colorbar.right + 5, y);
        RECT tick_label{colorbar.right + 8, y - 9, rect.right - 8, y + 10};
        DrawTextW(hdc, plot_number_wide(value).c_str(), -1, &tick_label, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
    }
    SelectObject(hdc, old_tick_pen); DeleteObject(tick_pen);
    SelectObject(hdc, g_font);

    RECT source_rect{rect.left + 16, plot.bottom + 50, rect.right - 16, plot.bottom + 70};
    SetTextColor(hdc, kColMuted);
    DrawTextW(hdc, g_preview.source_path.c_str(), -1, &source_rect, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
}

// ----- input model (sections + rows) ---------------------------------------
int g_section_geometry = 0, g_section_inlet = 0;

void add_edit_row(int section, int edit_id, int unit_id, const std::vector<const wchar_t*>& units,
                  const wchar_t* label, const wchar_t* suffix, bool advanced = false,
                  ParamVisibility vis = ParamVisibility::Always) {
    Row r;
    r.kind = RowKind::Edit;
    r.section = section;
    r.advanced = advanced;
    r.visibility = vis;
    r.value_id = edit_id;
    r.unit_id = unit_id;
    r.label = label;
    if (suffix) r.suffix = suffix;
    r.value = create_edit(g_rail, edit_id);
    if (unit_id) r.unit = create_combo(g_rail, unit_id, units);
    g_rows.push_back(r);
}

void add_combo_row(int section, int combo_id, const std::vector<const wchar_t*>& items,
                   const wchar_t* label, bool advanced = false) {
    Row r;
    r.kind = RowKind::Combo;
    r.section = section;
    r.advanced = advanced;
    r.value_id = combo_id;
    r.label = label;
    r.value = create_combo(g_rail, combo_id, items);
    g_rows.push_back(r);
}

void add_check_row(int section, int id, const wchar_t* text, bool advanced = false) {
    Row r;
    r.kind = RowKind::Check;
    r.section = section;
    r.advanced = advanced;
    r.value_id = id;
    r.value = create_checkbox(g_rail, id, text);
    g_rows.push_back(r);
}

int add_section(const wchar_t* title, bool collapsed = false, bool is_raw = false) {
    Section s;
    s.title = title;
    s.collapsed = collapsed;
    s.is_raw = is_raw;
    g_sections.push_back(s);
    return static_cast<int>(g_sections.size()) - 1;
}

void build_model() {
    const std::vector<const wchar_t*> len{L"m", L"mm", L"cm", L"um"};
    const std::vector<const wchar_t*> ang{L"deg", L"rad"};
    const std::vector<const wchar_t*> tilt{L"m/m", L"deg", L"rad"};
    const std::vector<const wchar_t*> tim{L"s", L"ms"};
    const std::vector<const wchar_t*> prs{L"Pa", L"kPa", L"MPa"};
    const std::vector<const wchar_t*> rpm{L"rad/s", L"rpm"};
    const std::vector<const wchar_t*> bc{L"DIRICHLET", L"NEUMANN", L"INLET_OUTLET"};
    const std::vector<const wchar_t*> methods{L"EULER_IMPLICIT", L"EULER_EXPLICIT", L"CRANK_NICOLSON", L"RK2", L"RK4"};

    g_section_geometry = add_section(L"Geometry");
    add_edit_row(g_section_geometry, kEditR, kComboUnitR, len, L"R", nullptr);
    add_edit_row(g_section_geometry, kEditC, kComboUnitC, len, L"c (clearance)", nullptr);
    add_edit_row(g_section_geometry, kEditE, kComboUnitE, len, L"e (eccentricity)", nullptr);
    add_edit_row(g_section_geometry, kEditL, kComboUnitL, len, L"L", nullptr);
    add_edit_row(g_section_geometry, kEditAttitude, kComboUnitAttitude, ang, L"attitude", nullptr, true);
    add_edit_row(g_section_geometry, kEditLoadAngle, kComboUnitLoadAngle, ang, L"load angle", nullptr, true);
    add_edit_row(g_section_geometry, kEditTiltX, kComboUnitTiltX, tilt, L"tilt_x", nullptr, true);
    add_edit_row(g_section_geometry, kEditTiltY, kComboUnitTiltY, tilt, L"tilt_y", nullptr, true);

    const int sec_op = add_section(L"Operating / Fluid");
    add_edit_row(sec_op, kEditOmega, kComboUnitOmega, rpm, L"omega", nullptr);
    add_edit_row(sec_op, kEditMu, 0, {}, L"mu", L"Pa s");
    add_edit_row(sec_op, kEditRho, 0, {}, L"rho", L"kg/m3");
    add_edit_row(sec_op, kEditPCav, kComboUnitPCav, prs, L"p_cav", nullptr, true);
    add_edit_row(sec_op, kEditBulkModulus, kComboUnitBulkModulus, prs, L"bulk modulus", nullptr, true);

    const int sec_energy = add_section(L"Energy");
    add_combo_row(sec_energy, kComboTemperatureModel, {L"ISOTHERMAL", L"ENERGY_EQUATION"}, L"thermal model");
    add_combo_row(sec_energy, kComboTemperatureTimeMethod, methods, L"thermal step", true);
    add_edit_row(sec_energy, kEditTemperatureInitial, 0, {}, L"T initial", L"K");
    add_edit_row(sec_energy, kEditJournalWallTemperature, 0, {}, L"T journal", L"K");
    add_edit_row(sec_energy, kEditBearingWallTemperature, 0, {}, L"T bearing", L"K");
    add_edit_row(sec_energy, kEditTemperatureReference, 0, {}, L"T reference", L"K", true);
    add_edit_row(sec_energy, kEditRhoCp, 0, {}, L"rho_cp", L"J/m3/K", true);
    add_edit_row(sec_energy, kEditThermalConductivity, 0, {}, L"k thermal", L"W/m/K", true);
    add_edit_row(sec_energy, kEditJournalHeatTransfer, 0, {}, L"h journal", L"W/m2/K", true);
    add_edit_row(sec_energy, kEditBearingHeatTransfer, 0, {}, L"h bearing", L"W/m2/K", true);

    const int sec_mesh = add_section(L"Mesh & Time");
    add_edit_row(sec_mesh, kEditNTheta, 0, {}, L"n_theta", L"cells");
    add_edit_row(sec_mesh, kEditNZ, 0, {}, L"n_z", L"cells");
    add_edit_row(sec_mesh, kEditEndT, kComboUnitEndT, tim, L"end_t", nullptr);
    add_edit_row(sec_mesh, kEditDt, kComboUnitDt, tim, L"dt", nullptr);
    add_edit_row(sec_mesh, kEditWriteInterval, kComboUnitWriteInterval, tim, L"write interval", nullptr);
    add_combo_row(sec_mesh, kComboSolutionMode, {L"TRANSIENT", L"STEADY_STATE"}, L"solution mode");

    const int sec_motion = add_section(L"Motion / Loads");
    add_combo_row(sec_motion, kComboMotionModel, {L"STATIC", L"MOVING_BEARING"}, L"motion model");
    add_combo_row(sec_motion, kComboPressureTimeMethod, methods, L"pressure step");
    add_combo_row(sec_motion, kComboMotionTimeMethod, methods, L"motion step");
    add_edit_row(sec_motion, kEditExternalLoadMagnitude, 0, {}, L"|F_ext,xy|", L"N");
    add_edit_row(sec_motion, kEditExternalLoadDirection, 0, {}, L"F_ext dir", L"deg");
    add_edit_row(sec_motion, kEditExternalLoadZ, 0, {}, L"F_ext z", L"N");
    add_check_row(sec_motion, kCheckBearingInitialFromAttitude, L"initialize from e and attitude", true);
    add_edit_row(sec_motion, kEditBearingInitialX, 0, {}, L"x0", L"m", true);
    add_edit_row(sec_motion, kEditBearingInitialY, 0, {}, L"y0", L"m", true);
    add_edit_row(sec_motion, kEditBearingInitialZ, 0, {}, L"z0", L"m", true);
    add_edit_row(sec_motion, kEditBearingInitialVx, 0, {}, L"vx0", L"m/s", true);
    add_edit_row(sec_motion, kEditBearingInitialVy, 0, {}, L"vy0", L"m/s", true);
    add_edit_row(sec_motion, kEditBearingInitialVz, 0, {}, L"vz0", L"m/s", true);
    add_edit_row(sec_motion, kEditBearingMass, 0, {}, L"mass", L"kg", true);
    add_edit_row(sec_motion, kEditBearingStiffnessX, 0, {}, L"k_x", L"N/m", true);
    add_edit_row(sec_motion, kEditBearingStiffnessY, 0, {}, L"k_y", L"N/m", true);
    add_edit_row(sec_motion, kEditBearingStiffnessZ, 0, {}, L"k_z", L"N/m", true);
    add_edit_row(sec_motion, kEditBearingDampingX, 0, {}, L"c_x", L"N s/m", true);
    add_edit_row(sec_motion, kEditBearingDampingY, 0, {}, L"c_y", L"N s/m", true);
    add_edit_row(sec_motion, kEditBearingDampingZ, 0, {}, L"c_z", L"N s/m", true);
    add_edit_row(sec_motion, kEditMinFilmThickness, 0, {}, L"h_min stop", L"m", true);
    add_check_row(sec_motion, kCheckStopOnNonpositiveFilm, L"stop on nonpositive film", true);

    const int sec_cav = add_section(L"Cavitation");
    add_combo_row(sec_cav, kComboCavitation, {L"ELROD_ADAMS", L"GUMBEL", L"FULL_SOMMERFELD"}, L"model");
    add_edit_row(sec_cav, kEditMaxOuter, 0, {}, L"max outer", L"iter", true);
    add_edit_row(sec_cav, kEditOuterTol, 0, {}, L"outer tol", L"-", true);
    add_edit_row(sec_cav, kEditThetaMin, 0, {}, L"theta min", L"-", true);
    add_check_row(sec_cav, kCheckLogOuter, L"log outer iterations", true);

    const int sec_bc = add_section(L"Axial Boundaries", true);
    add_combo_row(sec_bc, kComboBcSouth, bc, L"south type");
    add_edit_row(sec_bc, kEditBcSouthVal, kComboUnitBcSouthVal, prs, L"south p", nullptr);
    add_edit_row(sec_bc, kEditBcSouthTheta, 0, {}, L"south film", L"-", true);
    add_combo_row(sec_bc, kComboBcNorth, bc, L"north type");
    add_edit_row(sec_bc, kEditBcNorthVal, kComboUnitBcNorthVal, prs, L"north p", nullptr);
    add_edit_row(sec_bc, kEditBcNorthTheta, 0, {}, L"north film", L"-", true);

    g_section_inlet = add_section(L"Inlet");
    add_combo_row(g_section_inlet, kComboInletType, {L"GROOVE", L"CIRCULAR"}, L"type");
    add_edit_row(g_section_inlet, kEditInletTheta, kComboUnitInletTheta, ang, L"theta", nullptr);
    add_edit_row(g_section_inlet, kEditInletZ, kComboUnitInletZ, len, L"z", nullptr, false, ParamVisibility::CircularInletOnly);
    add_edit_row(g_section_inlet, kEditInletSize, kComboUnitInletSize, ang, L"width", nullptr);
    add_edit_row(g_section_inlet, kEditInletPressure, kComboUnitInletPressure, prs, L"p_supply", nullptr);

    const int sec_out = add_section(L"Output", true);
    add_edit_row(sec_out, kEditOutputDir, 0, {}, L"directory", nullptr);
    add_edit_row(sec_out, kEditFilenamePrefix, 0, {}, L"prefix", nullptr);
    add_check_row(sec_out, kCheckOutput3d, L"write curved 3D VTK/PVD", true);
    add_check_row(sec_out, kCheckOutputFlat, L"write flat unwrapped VTK/PVD", true);
    for (const auto& f : kOutputFields) add_check_row(sec_out, f.id, f.label, true);

    const int sec_raw = add_section(L"Raw config.txt", true, true);
    (void)sec_raw;
    create_button(g_rail, kSyncRawId, L"Sync from form");
    create_button(g_rail, kApplyRawId, L"Apply raw");
    HWND raw = create_child(g_rail, L"EDIT", L"",
        WS_TABSTOP | ES_MULTILINE | ES_AUTOVSCROLL | ES_AUTOHSCROLL | ES_WANTRETURN | WS_VSCROLL | WS_HSCROLL,
        WS_EX_CLIENTEDGE, kEditRawConfig, g_mono_font);
    (void)raw;
}

// ----- rail layout + paint -------------------------------------------------
bool row_advanced_hidden(const Row& r) { return r.advanced && !get_check(kCheckShowAdvanced); }

bool row_matches_filter(const Row& r, const std::wstring& filter_lower) {
    if (filter_lower.empty()) return true;
    std::wstring text = r.label;
    if (text.empty() && r.value != nullptr) text = get_window_text(r.value);
    std::transform(text.begin(), text.end(), text.begin(), ::towlower);
    return text.find(filter_lower) != std::wstring::npos;
}

bool row_visible(const Row& r, const std::wstring& filter_lower) {
    if (r.visibility == ParamVisibility::CircularInletOnly && !combo_is(kComboInletType, L"CIRCULAR")) return false;
    if (!filter_lower.empty()) return row_matches_filter(r, filter_lower);
    if (row_advanced_hidden(r)) return false;
    return true;
}

void show_row(Row& r, bool visible) {
    r.visible = visible;
    const int cmd = visible ? SW_SHOW : SW_HIDE;
    if (r.value) ShowWindow(r.value, cmd);
    if (r.unit) ShowWindow(r.unit, cmd);
}

void layout_rail() {
    if (g_rail == nullptr) return;
    RECT rc{};
    GetClientRect(g_rail, &rc);
    const int width = rc.right - rc.left;

    std::wstring filter = get_window_text(find_control(kEditSearch));
    std::transform(filter.begin(), filter.end(), filter.begin(), ::towlower);
    const bool searching = !filter.empty();

    const int pad = 12;
    const int label_x = pad + 4;
    const int label_w = 104;
    const int value_x = label_x + label_w + 6;
    const int unit_w = 62;
    const int gap = 6;
    const int value_w = std::max(70, width - value_x - unit_w - gap - pad);
    const int combo_w = value_w + gap + unit_w;
    const int row_h = 28;
    const int header_h = 30;

    int y = pad - g_rail_scroll_y;

    for (size_t si = 0; si < g_sections.size(); ++si) {
        Section& sec = g_sections[si];
        // Does this section have any visible row?
        bool any = false;
        for (Row& r : g_rows) if (r.section == static_cast<int>(si) && row_visible(r, filter)) { any = true; break; }
        if (sec.is_raw && (searching || get_check(kCheckShowAdvanced) || true)) {
            // Raw section header is always available (its body holds the editor).
            if (searching && !any) {
                // still show raw header so the editor stays reachable while searching? hide to reduce noise
            }
        }
        sec.visible = sec.is_raw ? !searching : any;
        if (!sec.visible) {
            // Hide all controls in a hidden section.
            for (Row& r : g_rows) if (r.section == static_cast<int>(si)) show_row(r, false);
            if (sec.is_raw) {
                ShowWindow(find_control(kSyncRawId), SW_HIDE);
                ShowWindow(find_control(kApplyRawId), SW_HIDE);
                ShowWindow(find_control(kEditRawConfig), SW_HIDE);
            }
            continue;
        }

        sec.header_rect = RECT{pad, y, width - pad, y + header_h};
        y += header_h;

        const bool expanded = searching || !sec.collapsed;

        if (!expanded) {
            for (Row& r : g_rows) if (r.section == static_cast<int>(si)) show_row(r, false);
            if (sec.is_raw) {
                ShowWindow(find_control(kSyncRawId), SW_HIDE);
                ShowWindow(find_control(kApplyRawId), SW_HIDE);
                ShowWindow(find_control(kEditRawConfig), SW_HIDE);
            }
            continue;
        }

        if (sec.is_raw) {
            move_control(find_control(kSyncRawId), label_x, y, 120, 26);
            move_control(find_control(kApplyRawId), label_x + 128, y, 100, 26);
            ShowWindow(find_control(kSyncRawId), SW_SHOW);
            ShowWindow(find_control(kApplyRawId), SW_SHOW);
            y += 32;
            const int raw_h = 320;
            move_control(find_control(kEditRawConfig), label_x, y, width - label_x - pad, raw_h);
            ShowWindow(find_control(kEditRawConfig), SW_SHOW);
            y += raw_h + 6;
            continue;
        }

        for (Row& r : g_rows) {
            if (r.section != static_cast<int>(si)) continue;
            const bool vis = row_visible(r, filter);
            show_row(r, vis);
            if (!vis) continue;
            if (r.kind == RowKind::Check) {
                move_control(r.value, label_x, y, width - label_x - pad, 22);
                r.label_rect = RECT{0, 0, 0, 0};
                y += row_h;
                continue;
            }
            r.label_rect = RECT{label_x, y + 4, value_x - 4, y + row_h - 2};
            if (r.kind == RowKind::Combo) {
                move_control(r.value, value_x, y - 1, combo_w, 240);
            } else {
                move_control(r.value, value_x, y, value_w, 24);
                if (r.unit) move_control(r.unit, value_x + value_w + gap, y - 1, unit_w, 240);
                else r.suffix_rect = RECT{value_x + value_w + gap, y + 4, value_x + combo_w, y + row_h - 2};
            }
            y += row_h;
        }
        y += 4;
    }

    g_rail_content_height = y + g_rail_scroll_y + pad;
    const int visible_h = rc.bottom - rc.top;
    const int max_scroll = std::max(0, g_rail_content_height - visible_h);
    g_rail_scroll_y = std::clamp(g_rail_scroll_y, 0, max_scroll);

    SCROLLINFO info{};
    info.cbSize = sizeof(info);
    info.fMask = SIF_RANGE | SIF_PAGE | SIF_POS;
    info.nMin = 0;
    info.nMax = std::max(g_rail_content_height, visible_h);
    info.nPage = static_cast<UINT>(std::max(1, visible_h));
    info.nPos = g_rail_scroll_y;
    SetScrollInfo(g_rail, SB_VERT, &info, TRUE);
    InvalidateRect(g_rail, nullptr, TRUE);
}

void scroll_rail(int delta) {
    if (g_rail == nullptr) return;
    RECT rc{};
    GetClientRect(g_rail, &rc);
    const int max_scroll = std::max(0, g_rail_content_height - static_cast<int>(rc.bottom - rc.top));
    const int next = std::clamp(g_rail_scroll_y + delta, 0, max_scroll);
    if (next == g_rail_scroll_y) return;
    g_rail_scroll_y = next;
    layout_rail();
}

void paint_rail(HDC hdc) {
    RECT rc{};
    GetClientRect(g_rail, &rc);
    FillRect(hdc, &rc, g_brush_page);
    SetBkMode(hdc, TRANSPARENT);

    std::wstring filter = get_window_text(find_control(kEditSearch));
    std::transform(filter.begin(), filter.end(), filter.begin(), ::towlower);
    const bool searching = !filter.empty();

    for (const Section& sec : g_sections) {
        if (!sec.visible) continue;
        // header band
        RECT hr = sec.header_rect;
        HPEN pen = CreatePen(PS_SOLID, 1, kColBorder);
        HGDIOBJ op = SelectObject(hdc, pen);
        MoveToEx(hdc, hr.left, hr.bottom - 1, nullptr);
        LineTo(hdc, hr.right, hr.bottom - 1);
        SelectObject(hdc, op);
        DeleteObject(pen);
        SelectObject(hdc, g_heading_font);
        SetTextColor(hdc, kColAccent);
        const bool expanded = searching || !sec.collapsed;
        const std::wstring marker = expanded ? L"\x25be  " : L"\x25b8  ";
        RECT tr{hr.left + 2, hr.top + 4, hr.right, hr.bottom};
        DrawTextW(hdc, (marker + sec.title).c_str(), -1, &tr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
    }

    SelectObject(hdc, g_font);
    for (const Row& r : g_rows) {
        if (!r.visible) continue;
        if (!r.label.empty() && r.label_rect.right > r.label_rect.left) {
            SetTextColor(hdc, kColText);
            RECT lr = r.label_rect;
            DrawTextW(hdc, r.label.c_str(), -1, &lr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
        }
        if (!r.suffix.empty() && r.suffix_rect.right > r.suffix_rect.left) {
            SetTextColor(hdc, kColMuted);
            RECT sr = r.suffix_rect;
            DrawTextW(hdc, r.suffix.c_str(), -1, &sr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
        }
    }
}

bool rail_header_hit(int x, int y) {
    for (Section& sec : g_sections) {
        if (!sec.visible) continue;
        if (x >= sec.header_rect.left && x < sec.header_rect.right &&
            y >= sec.header_rect.top && y < sec.header_rect.bottom) {
            sec.collapsed = !sec.collapsed;
            layout_rail();
            return true;
        }
    }
    return false;
}

// ----- summary card paint --------------------------------------------------
void paint_summary(HDC hdc) {
    RECT rc{};
    GetClientRect(g_summary, &rc);
    FillRect(hdc, &rc, g_brush_card);
    HPEN pen = CreatePen(PS_SOLID, 1, kColBorder);
    HGDIOBJ op = SelectObject(hdc, pen);
    HGDIOBJ ob = SelectObject(hdc, GetStockObject(NULL_BRUSH));
    Rectangle(hdc, rc.left, rc.top, rc.right, rc.bottom);
    SelectObject(hdc, ob);
    SelectObject(hdc, op);
    DeleteObject(pen);

    SetBkMode(hdc, TRANSPARENT);
    SelectObject(hdc, g_heading_font);
    SetTextColor(hdc, kColAccent);
    RECT tr{rc.left + 12, rc.top + 6, rc.right - 12, rc.top + 24};
    DrawTextW(hdc, L"LIVE SUMMARY", -1, &tr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);

    SelectObject(hdc, g_small_font);
    int y = rc.top + 28;
    const int col_w = (rc.right - rc.left - 24) / 2;
    int idx = 0;
    for (const auto& pair : g_summary_pairs) {
        const int col = idx % 2;
        const int row = idx / 2;
        const int x = rc.left + 12 + col * col_w;
        const int ry = y + row * 18;
        if (ry > rc.bottom - 38) break;
        SetTextColor(hdc, kColMuted);
        RECT lr{x, ry, x + 66, ry + 16};
        DrawTextW(hdc, pair.first.c_str(), -1, &lr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_NOPREFIX);
        SetTextColor(hdc, kColText);
        RECT vr{x + 66, ry, x + col_w - 6, ry + 16};
        DrawTextW(hdc, pair.second.c_str(), -1, &vr, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
        ++idx;
    }

    // Status line at the bottom.
    RECT sr{rc.left + 12, rc.bottom - 34, rc.right - 12, rc.bottom - 6};
    if (g_problems.empty()) {
        SetTextColor(hdc, kColOk);
        DrawTextW(hdc, L"\x2713 Ready to run.", -1, &sr, DT_LEFT | DT_TOP | DT_WORDBREAK | DT_NOPREFIX);
    } else {
        SetTextColor(hdc, g_has_errors ? kColInvalidText : kColWarn);
        std::wstring msg = g_problems.front();
        if (g_problems.size() > 1) msg += L"   (+" + std::to_wstring(g_problems.size() - 1) + L" more)";
        DrawTextW(hdc, msg.c_str(), -1, &sr, DT_LEFT | DT_TOP | DT_WORDBREAK | DT_NOPREFIX);
    }
}

// ----- fonts ---------------------------------------------------------------
void create_fonts() {
    g_font = CreateFontW(-15, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    g_heading_font = CreateFontW(-15, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    g_small_font = CreateFontW(-13, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    g_mono_font = CreateFontW(-14, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_MODERN, L"Consolas");
}

HWND g_ranks_label = nullptr, g_field_label = nullptr, g_step_label = nullptr, g_console_label = nullptr;

void create_controls() {
    // Action bar: run/stop (owner-drawn), ranks, progress, file buttons.
    g_run_button = CreateWindowExW(0, L"BUTTON", L"Run",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_TABSTOP | BS_OWNERDRAW,
        0, 0, 10, 10, g_hwnd, reinterpret_cast<HMENU>(static_cast<INT_PTR>(kRunId)), GetModuleHandleW(nullptr), nullptr);
    SendMessageW(g_run_button, WM_SETFONT, reinterpret_cast<WPARAM>(g_heading_font), TRUE);
    remember_control(kRunId, g_run_button);
    g_stop_button = CreateWindowExW(0, L"BUTTON", L"Stop",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_TABSTOP | BS_OWNERDRAW,
        0, 0, 10, 10, g_hwnd, reinterpret_cast<HMENU>(static_cast<INT_PTR>(kStopId)), GetModuleHandleW(nullptr), nullptr);
    SendMessageW(g_stop_button, WM_SETFONT, reinterpret_cast<WPARAM>(g_heading_font), TRUE);
    remember_control(kStopId, g_stop_button);
    EnableWindow(g_stop_button, FALSE);

    g_ranks_label = create_child(g_hwnd, L"STATIC", L"ranks", WS_VISIBLE | SS_LEFT, 0, 0);
    create_edit(g_hwnd, kEditMpiRanks);
    ShowWindow(find_control(kEditMpiRanks), SW_SHOW);

    g_progress = CreateWindowExW(0, PROGRESS_CLASSW, L"",
        WS_CHILD | WS_CLIPSIBLINGS, 0, 0, 10, 8, g_hwnd,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kProgressId)), GetModuleHandleW(nullptr), nullptr);
    SendMessageW(g_progress, PBM_SETRANGE, 0, MAKELPARAM(0, 100));

    create_button(g_hwnd, kOpenId, L"Open…");
    create_button(g_hwnd, kSaveId, L"Save");
    create_button(g_hwnd, kSaveAsId, L"Save As…");
    create_button(g_hwnd, kReloadId, L"Reload");
    create_button(g_hwnd, kDefaultsId, L"Reset");
    create_button(g_hwnd, kSolverBrowseId, L"Solver…");

    // Rail toolbar: search + advanced + preset + unit system (children of main
    // window so they stay fixed while the rail scrolls).
    HWND search = create_edit(g_hwnd, kEditSearch);
    ShowWindow(search, SW_SHOW);
    SendMessageW(search, EM_SETCUEBANNER, TRUE, reinterpret_cast<LPARAM>(L"Search parameters…"));
    HWND adv = create_checkbox(g_hwnd, kCheckShowAdvanced, L"Advanced");
    ShowWindow(adv, SW_SHOW);
    HWND preset = create_combo(g_hwnd, kComboPreset,
        {L"Load preset…", L"Default bearing", L"High eccentricity", L"Misaligned", L"Circular oil hole"});
    ShowWindow(preset, SW_SHOW);
    HWND units = create_combo(g_hwnd, kComboUnitSystem, {L"SI (m, Pa, rad/s)", L"Engineering (mm, MPa, rpm)"});
    ShowWindow(units, SW_SHOW);

    // Results toolbar + canvas.
    g_field_label = create_child(g_hwnd, L"STATIC", L"field", WS_VISIBLE | SS_LEFT, 0, 0);
    HWND fcombo = create_combo(g_hwnd, kComboPreviewField, {L"pressure"});
    ShowWindow(fcombo, SW_SHOW);
    g_step_label = create_child(g_hwnd, L"STATIC", L"step", WS_VISIBLE | SS_LEFT, 0, 0);
    HWND scombo = create_combo(g_hwnd, kComboPreviewStep, {L"Latest"});
    ShowWindow(scombo, SW_SHOW);
    create_button(g_hwnd, kPrevStepId, L"◀");
    create_button(g_hwnd, kNextStepId, L"▶");
    create_button(g_hwnd, kRefreshPreviewId, L"Refresh");
    create_button(g_hwnd, kOpenResultsId, L"Open folder");
    create_button(g_hwnd, kOpenParaviewId, L"Open in ParaView");

    g_preview_canvas = CreateWindowExW(WS_EX_CLIENTEDGE, kPreviewClass, L"",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS, 0, 0, 10, 10, g_hwnd,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kPreviewCanvasId)), GetModuleHandleW(nullptr), nullptr);

    // Console.
    g_console_label = create_child(g_hwnd, L"STATIC", L"CONSOLE", WS_VISIBLE | SS_LEFT, 0, 0, g_heading_font);
    create_button(g_hwnd, kClearLogId, L"Clear");
    HWND log = create_child(g_hwnd,
        g_richedit_available ? L"RICHEDIT50W" : L"EDIT", L"",
        WS_VISIBLE | WS_TABSTOP | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY | WS_VSCROLL,
        WS_EX_CLIENTEDGE, kEditLog, g_mono_font);
    if (g_richedit_available) SendMessageW(log, EM_SETBKGNDCOLOR, 0, RGB(250, 251, 253));

    // Rail panel (scrollable) must exist before the rows are created in it.
    g_rail = CreateWindowExW(WS_EX_CLIENTEDGE, kRailClass, L"",
        WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_CLIPCHILDREN, 0, 0, 10, 10, g_hwnd,
        nullptr, GetModuleHandleW(nullptr), nullptr);
    g_summary = CreateWindowExW(0, kSummaryClass, L"",
        WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS, 0, 0, 10, 10, g_hwnd,
        nullptr, GetModuleHandleW(nullptr), nullptr);

    g_rail_splitter = CreateWindowExW(0, kSplitterClass, L"", WS_CHILD | WS_VISIBLE,
        0, 0, 10, 10, g_hwnd, nullptr, GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(g_rail_splitter, GWLP_USERDATA, 0); // vertical
    g_console_splitter = CreateWindowExW(0, kSplitterClass, L"", WS_CHILD | WS_VISIBLE,
        0, 0, 10, 10, g_hwnd, nullptr, GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(g_console_splitter, GWLP_USERDATA, 1); // horizontal

    build_model();
    populate_preview_fields();
}

// ----- top-level layout ----------------------------------------------------
void layout_everything() {
    if (g_hwnd == nullptr) return;
    RECT client{};
    GetClientRect(g_hwnd, &client);
    const int W = client.right - client.left;
    const int H = client.bottom - client.top;

    // Action bar.
    move_control(g_run_button, kMargin, 11, 66, 34);
    move_control(g_stop_button, kMargin + 72, 11, 66, 34);
    move_control(g_ranks_label, kMargin + 150, 19, 38, 20);
    move_control_id(kEditMpiRanks, kMargin + 190, 14, 44, 26);
    move_control(g_progress, kMargin + 248, 40, 196, 8);

    int rx_btn = W - kMargin;
    auto place_right = [&](int id, int w) { rx_btn -= w; move_control_id(id, rx_btn, 12, w, 30); rx_btn -= 6; };
    place_right(kSolverBrowseId, 74);
    place_right(kDefaultsId, 64);
    place_right(kReloadId, 66);
    place_right(kSaveAsId, 84);
    place_right(kSaveId, 58);
    place_right(kOpenId, 66);

    const int top = kActionBarH + kGap;
    const int rail_x = kMargin;
    const int rail_w = g_rail_width;
    const int toolbar_h = 70;

    move_control_id(kEditSearch, rail_x, top, rail_w - 104, 26);
    move_control_id(kCheckShowAdvanced, rail_x + rail_w - 98, top + 2, 98, 22);
    const int half = (rail_w - 6) / 2;
    move_control_id(kComboPreset, rail_x, top + 34, half, 240);
    move_control_id(kComboUnitSystem, rail_x + half + 6, top + 34, half, 240);

    const int summary_y = H - kMargin - kSummaryH;
    move_control(g_summary, rail_x, summary_y, rail_w, kSummaryH);
    const int panel_y = top + toolbar_h + 4;
    move_control(g_rail, rail_x, panel_y, rail_w, std::max(60, summary_y - kGap - panel_y));

    move_control(g_rail_splitter, rail_x + rail_w, top, kSplitterW, H - top - kMargin);

    // Right column.
    const int rx = rail_x + rail_w + kSplitterW + kGap;
    const int rw = std::max(360, W - rx - kMargin);
    const int r_bottom = H - kMargin;
    g_console_height = std::clamp(g_console_height, kConsoleMinHeight, std::max(kConsoleMinHeight, (r_bottom - top) - kViewportMinHeight));
    const int csplit_y = r_bottom - g_console_height - kSplitterW;
    move_control(g_console_splitter, rx, csplit_y, rw, kSplitterW);

    // Results toolbar.
    move_control(g_field_label, rx, top + 6, 32, 20);
    move_control_id(kComboPreviewField, rx + 36, top + 2, 150, 240);
    move_control(g_step_label, rx + 194, top + 6, 28, 20);
    move_control_id(kComboPreviewStep, rx + 224, top + 2, 108, 240);
    move_control_id(kPrevStepId, rx + 338, top + 2, 30, 26);
    move_control_id(kNextStepId, rx + 370, top + 2, 30, 26);
    move_control_id(kRefreshPreviewId, rx + 408, top + 2, 78, 26);
    move_control_id(kOpenResultsId, rx + 492, top + 2, 96, 26);
    move_control_id(kOpenParaviewId, rx + 594, top + 2, 120, 26);

    const int canvas_top = top + 36;
    move_control(g_preview_canvas, rx, canvas_top, rw, std::max(60, csplit_y - kGap - canvas_top));

    // Console.
    const int console_top = csplit_y + kSplitterW + kGap;
    move_control(g_console_label, rx, console_top, 100, 18);
    move_control_id(kClearLogId, rx + rw - 64, console_top - 2, 64, 24);
    move_control_id(kEditLog, rx, console_top + 22, rw, std::max(40, r_bottom - (console_top + 22)));

    layout_rail();
    InvalidateRect(g_hwnd, nullptr, FALSE);
}

void update_rail_splitter(int client_x) {
    RECT client{};
    GetClientRect(g_hwnd, &client);
    const int max_w = std::max(kRailMinWidth, static_cast<int>(client.right - client.left) - kRailMaxSlack);
    g_rail_width = std::clamp(client_x - kMargin, kRailMinWidth, max_w);
    layout_everything();
}

void update_console_splitter(int client_y) {
    RECT client{};
    GetClientRect(g_hwnd, &client);
    const int H = client.bottom - client.top;
    const int top = kActionBarH + kGap;
    const int r_bottom = H - kMargin;
    const int requested = r_bottom - client_y - kSplitterW;
    g_console_height = std::clamp(requested, kConsoleMinHeight, std::max(kConsoleMinHeight, (r_bottom - top) - kViewportMinHeight));
    layout_everything();
}

// ----- action-bar paint + owner-drawn buttons ------------------------------
void draw_action_button(LPDRAWITEMSTRUCT d, COLORREF base, COLORREF hot, const wchar_t* text) {
    RECT r = d->rcItem;
    const bool disabled = (d->itemState & ODS_DISABLED) != 0;
    const bool pressed = (d->itemState & ODS_SELECTED) != 0;
    const COLORREF c = disabled ? RGB(201, 207, 215) : (pressed ? hot : base);
    HBRUSH b = CreateSolidBrush(c);
    FillRect(d->hDC, &r, b);
    DeleteObject(b);
    SetBkMode(d->hDC, TRANSPARENT);
    SelectObject(d->hDC, g_heading_font);
    SetTextColor(d->hDC, disabled ? RGB(244, 247, 250) : RGB(255, 255, 255));
    DrawTextW(d->hDC, text, -1, &r, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
}

void paint_action_bar(HDC hdc, const RECT& client) {
    HPEN pen = CreatePen(PS_SOLID, 1, kColBorder);
    HGDIOBJ op = SelectObject(hdc, pen);
    MoveToEx(hdc, client.left, kActionBarH - 1, nullptr);
    LineTo(hdc, client.right, kActionBarH - 1);
    SelectObject(hdc, op);
    DeleteObject(pen);

    // Run-state chip.
    RECT chip{kMargin + 248, 12, kMargin + 248 + 196, 12 + 24};
    if (g_run_state == RunState::Running) chip.bottom = 12 + 22; // leave room for the progress bar below
    COLORREF cc = kColMuted;
    std::wstring text = L"Idle";
    std::ostringstream o;
    switch (g_run_state) {
        case RunState::Running:
            cc = kColAccent;
            o << "Running  " << static_cast<int>(g_run_progress * 100.0) << "%  \xC2\xB7  " << std::fixed << std::setprecision(1) << g_run_elapsed_s << "s";
            text = utf8_to_wide(o.str());
            break;
        case RunState::Done:
            cc = kColOk;
            o << "Done  \xC2\xB7  " << std::fixed << std::setprecision(1) << g_run_elapsed_s << "s";
            text = utf8_to_wide(o.str());
            break;
        case RunState::Failed: cc = kColStop; text = L"Failed"; break;
        default: break;
    }
    HBRUSH cb = CreateSolidBrush(cc);
    FillRect(hdc, &chip, cb);
    DeleteObject(cb);
    SetBkMode(hdc, TRANSPARENT);
    SelectObject(hdc, g_small_font);
    SetTextColor(hdc, RGB(255, 255, 255));
    RECT ct = chip; ct.left += 10;
    DrawTextW(hdc, text.c_str(), -1, &ct, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
}

// ----- shared control-colour handling --------------------------------------
LRESULT handle_ctl_color(UINT message, HDC dc, HWND control) {
    if (message == WM_CTLCOLOREDIT) {
        const int id = GetDlgCtrlID(control);
        if (g_invalid_ids.count(id)) {
            SetBkColor(dc, kColInvalidBg);
            SetTextColor(dc, kColInvalidText);
            return reinterpret_cast<LRESULT>(g_brush_invalid);
        }
        SetBkColor(dc, kColField);
        SetTextColor(dc, kColText);
        return reinterpret_cast<LRESULT>(g_brush_field);
    }
    if (message == WM_CTLCOLORLISTBOX) {
        SetBkColor(dc, kColField);
        SetTextColor(dc, kColText);
        return reinterpret_cast<LRESULT>(g_brush_field);
    }
    // WM_CTLCOLORSTATIC / BTN
    SetBkColor(dc, kColPage);
    SetTextColor(dc, kColText);
    return reinterpret_cast<LRESULT>(g_brush_page);
}

void on_input_changed(int id, int code) {
    if (g_loading) return;
    if (id == kEditRawConfig) {
        if (code == EN_CHANGE) { g_active_raw = true; mark_dirty(); }
        return;
    }
    if (id == kComboInletType && code == CBN_SELCHANGE) {
        update_inlet_editor_visibility();
        layout_rail();
    }
    g_active_raw = false;
    mark_dirty();
    recompute_summary();
}

// ----- window procedures ---------------------------------------------------
LRESULT CALLBACK rail_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_ERASEBKGND: {
            HDC hdc = reinterpret_cast<HDC>(wparam);
            RECT rc{};
            GetClientRect(hwnd, &rc);
            FillRect(hdc, &rc, g_brush_page);
            return 1;
        }
        case WM_PAINT: {
            PAINTSTRUCT ps{};
            HDC hdc = BeginPaint(hwnd, &ps);
            paint_rail(hdc);
            EndPaint(hwnd, &ps);
            return 0;
        }
        case WM_CTLCOLOREDIT:
        case WM_CTLCOLORLISTBOX:
        case WM_CTLCOLORSTATIC:
        case WM_CTLCOLORBTN:
            return handle_ctl_color(message, reinterpret_cast<HDC>(wparam), reinterpret_cast<HWND>(lparam));
        case WM_LBUTTONDOWN:
            SetFocus(hwnd);
            rail_header_hit(GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam));
            return 0;
        case WM_MOUSEWHEEL:
            scroll_rail(-GET_WHEEL_DELTA_WPARAM(wparam) / WHEEL_DELTA * 60);
            return 0;
        case WM_VSCROLL: {
            RECT rc{};
            GetClientRect(hwnd, &rc);
            const int page = std::max(1, static_cast<int>(rc.bottom - rc.top));
            switch (LOWORD(wparam)) {
                case SB_LINEUP: scroll_rail(-28); break;
                case SB_LINEDOWN: scroll_rail(28); break;
                case SB_PAGEUP: scroll_rail(-page); break;
                case SB_PAGEDOWN: scroll_rail(page); break;
                case SB_THUMBTRACK:
                case SB_THUMBPOSITION: {
                    SCROLLINFO info{};
                    info.cbSize = sizeof(info);
                    info.fMask = SIF_TRACKPOS;
                    GetScrollInfo(hwnd, SB_VERT, &info);
                    scroll_rail(info.nTrackPos - g_rail_scroll_y);
                    break;
                }
            }
            return 0;
        }
        case WM_COMMAND: {
            const int id = LOWORD(wparam);
            const int code = HIWORD(wparam);
            if (id == kSyncRawId && code == BN_CLICKED) { sync_raw_from_form(); return 0; }
            if (id == kApplyRawId && code == BN_CLICKED) { apply_raw_to_form(); return 0; }
            if (id == kEditRawConfig && code == EN_CHANGE) { on_input_changed(id, code); return 0; }
            if ((code == EN_CHANGE || code == CBN_SELCHANGE || code == BN_CLICKED) &&
                id >= kEditR && id <= kCheckFieldTorque) {
                on_input_changed(id, code);
                return 0;
            }
            break;
        }
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK summary_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    if (message == WM_PAINT) {
        PAINTSTRUCT ps{};
        HDC hdc = BeginPaint(hwnd, &ps);
        paint_summary(hdc);
        EndPaint(hwnd, &ps);
        return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK preview_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    if (message == WM_PAINT) {
        PAINTSTRUCT ps{};
        HDC hdc = BeginPaint(hwnd, &ps);
        draw_preview(hwnd, hdc);
        EndPaint(hwnd, &ps);
        return 0;
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK splitter_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    const bool horizontal = GetWindowLongPtrW(hwnd, GWLP_USERDATA) == 1;
    bool& dragging = horizontal ? g_console_splitter_dragging : g_rail_splitter_dragging;
    switch (message) {
        case WM_SETCURSOR:
            SetCursor(LoadCursorW(nullptr, horizontal ? IDC_SIZENS : IDC_SIZEWE));
            return TRUE;
        case WM_LBUTTONDOWN:
            dragging = true;
            SetCapture(hwnd);
            return 0;
        case WM_MOUSEMOVE:
            if (dragging) {
                POINT pt{};
                GetCursorPos(&pt);
                ScreenToClient(g_hwnd, &pt);
                if (horizontal) update_console_splitter(pt.y);
                else update_rail_splitter(pt.x);
            }
            return 0;
        case WM_LBUTTONUP:
            if (dragging) { dragging = false; ReleaseCapture(); }
            return 0;
        case WM_CAPTURECHANGED:
            dragging = false;
            return 0;
        case WM_PAINT: {
            PAINTSTRUCT ps{};
            HDC hdc = BeginPaint(hwnd, &ps);
            RECT rc{};
            GetClientRect(hwnd, &rc);
            FillRect(hdc, &rc, g_brush_page);
            HPEN pen = CreatePen(PS_SOLID, 1, kColBorder);
            HGDIOBJ op = SelectObject(hdc, pen);
            if (horizontal) {
                const int cy = (rc.bottom - rc.top) / 2;
                MoveToEx(hdc, rc.left + 16, cy, nullptr);
                LineTo(hdc, rc.right - 16, cy);
            } else {
                const int cx = (rc.right - rc.left) / 2;
                MoveToEx(hdc, cx, rc.top + 16, nullptr);
                LineTo(hdc, cx, rc.bottom - 16);
            }
            SelectObject(hdc, op);
            DeleteObject(pen);
            EndPaint(hwnd, &ps);
            return 0;
        }
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

LRESULT CALLBACK window_proc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam) {
    switch (message) {
        case WM_CREATE:
            g_hwnd = hwnd;
            g_executable_dir = executable_directory();
            g_config_path = g_executable_dir / L"config.txt";
            g_solver_path = (g_executable_dir / L"pancake.exe").wstring();
            create_fonts();
            create_controls();
            load_config();
            layout_everything();
            update_title();
            recompute_summary();
            return 0;

        case WM_SIZE:
            layout_everything();
            return 0;

        case WM_GETMINMAXINFO: {
            auto* mm = reinterpret_cast<MINMAXINFO*>(lparam);
            mm->ptMinTrackSize.x = 1120;
            mm->ptMinTrackSize.y = 700;
            return 0;
        }

        case WM_ERASEBKGND: {
            HDC hdc = reinterpret_cast<HDC>(wparam);
            RECT rc{};
            GetClientRect(hwnd, &rc);
            FillRect(hdc, &rc, g_brush_page);
            return 1;
        }

        case WM_PAINT: {
            PAINTSTRUCT ps{};
            HDC hdc = BeginPaint(hwnd, &ps);
            RECT rc{};
            GetClientRect(hwnd, &rc);
            paint_action_bar(hdc, rc);
            EndPaint(hwnd, &ps);
            return 0;
        }

        case WM_CTLCOLOREDIT:
        case WM_CTLCOLORLISTBOX:
        case WM_CTLCOLORSTATIC:
        case WM_CTLCOLORBTN:
            return handle_ctl_color(message, reinterpret_cast<HDC>(wparam), reinterpret_cast<HWND>(lparam));

        case WM_DRAWITEM: {
            auto* d = reinterpret_cast<LPDRAWITEMSTRUCT>(lparam);
            if (d->CtlID == static_cast<UINT>(kRunId)) { draw_action_button(d, kColRun, kColRunHot, L"▶  Run"); return TRUE; }
            if (d->CtlID == static_cast<UINT>(kStopId)) { draw_action_button(d, kColStop, kColStopHot, L"■  Stop"); return TRUE; }
            break;
        }

        case WM_TIMER:
            if (wparam == kProcessTimer) { poll_solver(); return 0; }
            break;

        case WM_COMMAND: {
            const int id = LOWORD(wparam);
            const int code = HIWORD(wparam);
            switch (id) {
                case kRunId: start_solver(); return 0;
                case kStopId: stop_solver(); return 0;
                case kSaveId: save_config(); return 0;
                case kSaveAsId: save_as_dialog(); return 0;
                case kOpenId: open_config_dialog(); return 0;
                case kReloadId: if (confirm_discard_changes()) load_config(); return 0;
                case kDefaultsId: reset_defaults(); return 0;
                case kSolverBrowseId: browse_solver(); return 0;
                case kRefreshPreviewId: if (save_config()) refresh_preview(true); return 0;
                case kOpenResultsId: open_results_folder(); return 0;
                case kOpenParaviewId: open_paraview_result(); return 0;
                case kPrevStepId: move_preview_step_selection(-1); return 0;
                case kNextStepId: move_preview_step_selection(1); return 0;
                case kClearLogId: clear_log(); return 0;
            }
            if (id == kComboPreset && code == CBN_SELCHANGE) {
                const int sel = static_cast<int>(SendMessageW(find_control(kComboPreset), CB_GETCURSEL, 0, 0));
                if (sel > 0) { load_preset(sel - 1); SendMessageW(find_control(kComboPreset), CB_SETCURSEL, 0, 0); }
                return 0;
            }
            if (id == kComboUnitSystem && code == CBN_SELCHANGE) {
                apply_unit_system(SendMessageW(find_control(kComboUnitSystem), CB_GETCURSEL, 0, 0) == 1);
                return 0;
            }
            if (id == kEditSearch && code == EN_CHANGE) { g_rail_scroll_y = 0; layout_rail(); return 0; }
            if (id == kCheckShowAdvanced && code == BN_CLICKED) { layout_rail(); return 0; }
            if ((id == kComboPreviewField || id == kComboPreviewStep) && code == CBN_SELCHANGE && !g_loading) {
                refresh_preview(false);
                return 0;
            }
            break;
        }

        case WM_CLOSE:
            if (g_solver.running) {
                const int result = MessageBoxW(hwnd, L"The solver is still running. Stop it and close?", L"Pancake",
                                               MB_ICONWARNING | MB_YESNO | MB_DEFBUTTON2);
                if (result != IDYES) return 0;
                stop_solver();
            }
            if (confirm_discard_changes()) DestroyWindow(hwnd);
            return 0;

        case WM_DESTROY:
            KillTimer(hwnd, kProcessTimer);
            close_solver_handles();
            if (g_font) DeleteObject(g_font);
            if (g_heading_font) DeleteObject(g_heading_font);
            if (g_small_font) DeleteObject(g_small_font);
            if (g_mono_font) DeleteObject(g_mono_font);
            if (g_brush_page) DeleteObject(g_brush_page);
            if (g_brush_field) DeleteObject(g_brush_field);
            if (g_brush_card) DeleteObject(g_brush_card);
            if (g_brush_invalid) DeleteObject(g_brush_invalid);
            if (g_richedit_module) FreeLibrary(g_richedit_module);
            PostQuitMessage(0);
            return 0;
    }
    return DefWindowProcW(hwnd, message, wparam, lparam);
}

} // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show_command) {
    INITCOMMONCONTROLSEX common{sizeof(INITCOMMONCONTROLSEX), ICC_STANDARD_CLASSES | ICC_PROGRESS_CLASS | ICC_TAB_CLASSES};
    InitCommonControlsEx(&common);
    g_richedit_module = LoadLibraryW(L"Msftedit.dll");
    g_richedit_available = g_richedit_module != nullptr;

    g_brush_page = CreateSolidBrush(kColPage);
    g_brush_field = CreateSolidBrush(kColField);
    g_brush_card = CreateSolidBrush(kColCard);
    g_brush_invalid = CreateSolidBrush(kColInvalidBg);

    auto register_simple = [&](const wchar_t* name, WNDPROC proc, HBRUSH bg, LPCWSTR cursor) {
        WNDCLASSEXW wc{};
        wc.cbSize = sizeof(wc);
        wc.lpfnWndProc = proc;
        wc.hInstance = instance;
        wc.hCursor = LoadCursorW(nullptr, cursor);
        wc.hbrBackground = bg;
        wc.lpszClassName = name;
        RegisterClassExW(&wc);
    };
    register_simple(kPreviewClass, preview_proc, g_brush_card, IDC_ARROW);
    register_simple(kSplitterClass, splitter_proc, g_brush_page, IDC_ARROW);
    register_simple(kRailClass, rail_proc, g_brush_page, IDC_ARROW);
    register_simple(kSummaryClass, summary_proc, g_brush_card, IDC_ARROW);

    HICON app_icon = static_cast<HICON>(LoadImageW(instance, MAKEINTRESOURCEW(IDI_PANCAKE), IMAGE_ICON, 0, 0, LR_DEFAULTSIZE | LR_SHARED));
    HICON app_icon_small = static_cast<HICON>(LoadImageW(instance, MAKEINTRESOURCEW(IDI_PANCAKE), IMAGE_ICON,
        GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), LR_SHARED));

    WNDCLASSEXW window_class{};
    window_class.cbSize = sizeof(window_class);
    window_class.lpfnWndProc = window_proc;
    window_class.hInstance = instance;
    window_class.hIcon = app_icon ? app_icon : LoadIconW(nullptr, IDI_APPLICATION);
    window_class.hIconSm = app_icon_small ? app_icon_small : window_class.hIcon;
    window_class.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    window_class.hbrBackground = g_brush_page;
    window_class.lpszClassName = kWindowClass;
    if (!RegisterClassExW(&window_class)) return 1;

    HWND hwnd = CreateWindowExW(0, kWindowClass, L"Pancake", WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 1280, 840, nullptr, nullptr, instance, nullptr);
    if (hwnd == nullptr) return 1;
    if (app_icon) SendMessageW(hwnd, WM_SETICON, ICON_BIG, reinterpret_cast<LPARAM>(app_icon));
    if (app_icon_small) SendMessageW(hwnd, WM_SETICON, ICON_SMALL, reinterpret_cast<LPARAM>(app_icon_small));

    ShowWindow(hwnd, show_command);
    UpdateWindow(hwnd);

    MSG message{};
    while (GetMessageW(&message, nullptr, 0, 0) > 0) {
        TranslateMessage(&message);
        DispatchMessageW(&message);
    }
    return static_cast<int>(message.wParam);
}
