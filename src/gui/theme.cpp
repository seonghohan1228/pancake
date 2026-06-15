#include "theme.hpp"

#include "embedded_jetbrains_mono.h"
#include "embedded_roboto_medium.h"
#include "embedded_roboto_regular.h"

namespace theme {

void apply_style(float dpi_scale) {
    ImGuiStyle style;  // fresh defaults so repeated rescaling never compounds
    ImGui::StyleColorsLight(&style);

    style.WindowRounding = 4.0f;
    style.ChildRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.PopupRounding = 4.0f;
    style.GrabRounding = 3.0f;
    style.TabRounding = 3.0f;
    style.ScrollbarRounding = 6.0f;
    style.WindowPadding = ImVec2(10, 10);
    style.FramePadding = ImVec2(8, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.ItemInnerSpacing = ImVec2(6, 4);
    style.IndentSpacing = 18.0f;
    style.CellPadding = ImVec2(6, 3);
    style.WindowBorderSize = 1.0f;
    style.FrameBorderSize = 1.0f;  // inputs need an outline on light backgrounds
    style.DockingSeparatorSize = 2.0f;

    ImVec4* colors = style.Colors;
    const ImVec4 bg0{0.953f, 0.961f, 0.973f, 1.0f};   // window background
    const ImVec4 bg1{0.910f, 0.922f, 0.937f, 1.0f};   // raised surfaces
    const ImVec4 bg2{0.855f, 0.871f, 0.892f, 1.0f};   // hovered
    const ImVec4 bg3{0.792f, 0.816f, 0.847f, 1.0f};   // active
    const ImVec4 field{1.0f, 1.0f, 1.0f, 1.0f};       // input fields
    const ImVec4 text{0.118f, 0.141f, 0.176f, 1.0f};

    colors[ImGuiCol_Text] = text;
    colors[ImGuiCol_TextDisabled] = kMuted;
    colors[ImGuiCol_WindowBg] = bg0;
    colors[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    colors[ImGuiCol_PopupBg] = ImVec4(1.0f, 1.0f, 1.0f, 0.99f);
    colors[ImGuiCol_Border] = ImVec4(0.69f, 0.72f, 0.76f, 0.60f);
    colors[ImGuiCol_FrameBg] = field;
    colors[ImGuiCol_FrameBgHovered] = ImVec4(0.95f, 0.97f, 1.0f, 1.0f);
    colors[ImGuiCol_FrameBgActive] = ImVec4(0.90f, 0.94f, 1.0f, 1.0f);
    colors[ImGuiCol_TitleBg] = bg1;
    colors[ImGuiCol_TitleBgActive] = bg2;
    colors[ImGuiCol_TitleBgCollapsed] = bg1;
    colors[ImGuiCol_MenuBarBg] = bg1;
    colors[ImGuiCol_ScrollbarBg] = bg0;
    colors[ImGuiCol_ScrollbarGrab] = bg2;
    colors[ImGuiCol_ScrollbarGrabHovered] = bg3;
    colors[ImGuiCol_ScrollbarGrabActive] = kAccentDim;
    colors[ImGuiCol_CheckMark] = kAccent;
    colors[ImGuiCol_SliderGrab] = kAccentDim;
    colors[ImGuiCol_SliderGrabActive] = kAccent;
    colors[ImGuiCol_Button] = bg1;
    colors[ImGuiCol_ButtonHovered] = bg2;
    colors[ImGuiCol_ButtonActive] = bg3;
    colors[ImGuiCol_Header] = bg1;
    colors[ImGuiCol_HeaderHovered] = bg2;
    colors[ImGuiCol_HeaderActive] = bg3;
    colors[ImGuiCol_Separator] = colors[ImGuiCol_Border];
    colors[ImGuiCol_SeparatorHovered] = kAccentDim;
    colors[ImGuiCol_SeparatorActive] = kAccent;
    colors[ImGuiCol_ResizeGrip] = bg2;
    colors[ImGuiCol_ResizeGripHovered] = bg3;
    colors[ImGuiCol_ResizeGripActive] = kAccent;
    colors[ImGuiCol_Tab] = bg1;
    colors[ImGuiCol_TabHovered] = bg2;
    colors[ImGuiCol_TabSelected] = bg0;
    colors[ImGuiCol_TabSelectedOverline] = kAccent;
    colors[ImGuiCol_TabDimmed] = bg1;
    colors[ImGuiCol_TabDimmedSelected] = bg0;
    colors[ImGuiCol_DockingPreview] = ImVec4(kAccent.x, kAccent.y, kAccent.z, 0.4f);
    colors[ImGuiCol_DockingEmptyBg] = bg0;
    colors[ImGuiCol_PlotLines] = kAccent;
    colors[ImGuiCol_PlotHistogram] = kAccent;
    colors[ImGuiCol_TableHeaderBg] = bg1;
    colors[ImGuiCol_TableBorderStrong] = colors[ImGuiCol_Border];
    colors[ImGuiCol_TableBorderLight] = ImVec4(0.78f, 0.80f, 0.83f, 0.5f);
    colors[ImGuiCol_TableRowBgAlt] = ImVec4(0, 0, 0, 0.025f);
    colors[ImGuiCol_TextSelectedBg] = ImVec4(kAccent.x, kAccent.y, kAccent.z, 0.30f);
    colors[ImGuiCol_NavCursor] = kAccent;
    colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0.2f, 0.2f, 0.2f, 0.35f);

    style.ScaleAllSizes(dpi_scale);
    ImGui::GetStyle() = style;
}

Fonts build_fonts(float dpi_scale) {
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();

    ImFontConfig config;
    config.FontDataOwnedByAtlas = false;  // data lives in the exe image

    Fonts fonts;
    fonts.sans = io.Fonts->AddFontFromMemoryTTF(
        const_cast<unsigned char*>(g_font_roboto_regular),
        static_cast<int>(g_font_roboto_regular_size), 16.0f * dpi_scale, &config);
    fonts.heading = io.Fonts->AddFontFromMemoryTTF(
        const_cast<unsigned char*>(g_font_roboto_medium),
        static_cast<int>(g_font_roboto_medium_size), 17.0f * dpi_scale, &config);
    fonts.mono = io.Fonts->AddFontFromMemoryTTF(
        const_cast<unsigned char*>(g_font_jetbrains_mono),
        static_cast<int>(g_font_jetbrains_mono_size), 14.5f * dpi_scale, &config);
    io.FontDefault = fonts.sans;
    io.Fonts->Build();
    return fonts;
}

}  // namespace theme
