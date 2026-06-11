#include "theme.hpp"

#include "embedded_jetbrains_mono.h"
#include "embedded_roboto_medium.h"
#include "embedded_roboto_regular.h"

namespace theme {

void apply_style(float dpi_scale) {
    ImGuiStyle style;  // fresh defaults so repeated rescaling never compounds
    ImGui::StyleColorsDark(&style);

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
    style.FrameBorderSize = 0.0f;
    style.DockingSeparatorSize = 2.0f;

    ImVec4* colors = style.Colors;
    const ImVec4 bg0{0.106f, 0.113f, 0.125f, 1.0f};   // window background
    const ImVec4 bg1{0.137f, 0.145f, 0.161f, 1.0f};   // child/frame background
    const ImVec4 bg2{0.180f, 0.190f, 0.208f, 1.0f};   // hovered
    const ImVec4 bg3{0.227f, 0.239f, 0.259f, 1.0f};   // active
    const ImVec4 text{0.886f, 0.898f, 0.914f, 1.0f};

    colors[ImGuiCol_Text] = text;
    colors[ImGuiCol_TextDisabled] = kMuted;
    colors[ImGuiCol_WindowBg] = bg0;
    colors[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    colors[ImGuiCol_PopupBg] = ImVec4(0.09f, 0.095f, 0.105f, 0.98f);
    colors[ImGuiCol_Border] = ImVec4(0.26f, 0.27f, 0.30f, 0.50f);
    colors[ImGuiCol_FrameBg] = bg1;
    colors[ImGuiCol_FrameBgHovered] = bg2;
    colors[ImGuiCol_FrameBgActive] = bg3;
    colors[ImGuiCol_TitleBg] = bg0;
    colors[ImGuiCol_TitleBgActive] = bg1;
    colors[ImGuiCol_TitleBgCollapsed] = bg0;
    colors[ImGuiCol_MenuBarBg] = bg1;
    colors[ImGuiCol_ScrollbarBg] = bg0;
    colors[ImGuiCol_ScrollbarGrab] = bg2;
    colors[ImGuiCol_ScrollbarGrabHovered] = bg3;
    colors[ImGuiCol_ScrollbarGrabActive] = kAccentDim;
    colors[ImGuiCol_CheckMark] = kAccent;
    colors[ImGuiCol_SliderGrab] = kAccentDim;
    colors[ImGuiCol_SliderGrabActive] = kAccent;
    colors[ImGuiCol_Button] = bg2;
    colors[ImGuiCol_ButtonHovered] = bg3;
    colors[ImGuiCol_ButtonActive] = kAccentDim;
    colors[ImGuiCol_Header] = bg2;
    colors[ImGuiCol_HeaderHovered] = bg3;
    colors[ImGuiCol_HeaderActive] = kAccentDim;
    colors[ImGuiCol_Separator] = colors[ImGuiCol_Border];
    colors[ImGuiCol_SeparatorHovered] = kAccentDim;
    colors[ImGuiCol_SeparatorActive] = kAccent;
    colors[ImGuiCol_ResizeGrip] = bg2;
    colors[ImGuiCol_ResizeGripHovered] = bg3;
    colors[ImGuiCol_ResizeGripActive] = kAccent;
    colors[ImGuiCol_Tab] = bg1;
    colors[ImGuiCol_TabHovered] = bg3;
    colors[ImGuiCol_TabSelected] = bg3;
    colors[ImGuiCol_TabSelectedOverline] = kAccent;
    colors[ImGuiCol_TabDimmed] = bg0;
    colors[ImGuiCol_TabDimmedSelected] = bg1;
    colors[ImGuiCol_DockingPreview] = ImVec4(kAccent.x, kAccent.y, kAccent.z, 0.5f);
    colors[ImGuiCol_DockingEmptyBg] = bg0;
    colors[ImGuiCol_PlotLines] = kAccent;
    colors[ImGuiCol_PlotHistogram] = kAccent;
    colors[ImGuiCol_TableHeaderBg] = bg1;
    colors[ImGuiCol_TableBorderStrong] = colors[ImGuiCol_Border];
    colors[ImGuiCol_TableBorderLight] = ImVec4(0.22f, 0.23f, 0.25f, 0.4f);
    colors[ImGuiCol_TableRowBgAlt] = ImVec4(1, 1, 1, 0.02f);
    colors[ImGuiCol_TextSelectedBg] = ImVec4(kAccent.x, kAccent.y, kAccent.z, 0.35f);
    colors[ImGuiCol_NavCursor] = kAccent;
    colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0, 0, 0, 0.5f);

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
