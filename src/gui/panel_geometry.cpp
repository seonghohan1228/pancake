// Geometry panel: bearing size, clearance, oil feed, and (Advanced) tilt.

#include "gui_app.hpp"
#include "imgui.h"

#include <cstdio>

/// Bearing diameter, length, radial clearance, oil feed; Advanced: hole position, tilt.
void GuiApp::draw_geometry_panel() {
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Geometry)];
    if (!ui::begin_form("##geometry")) return;
    bool changed = false;

    changed |= ui::form_number("Bearing diameter D", "Nominal bore (shaft) diameter, D = 2R. Typical: 15 to 40 mm.",
                               "geom.diameter", case_.bearing_diameter, gui::UnitFamily::Length, units());
    changed |= ui::form_number("Bearing length L", "Axial length of the film. Typical L/D: 0.5 to 1.2.",
                               "geom.length", case_.bearing_length, gui::UnitFamily::Length, units());
    changed |= ui::form_number("Radial clearance c",
                               "Bore radius minus shaft radius. Typical c/R: 0.001 to 0.002 (10 to 40 um).",
                               "geom.clearance", case_.radial_clearance, gui::UnitFamily::SmallLength,
                               units());

    changed |= ui::form_radio<gui::FeedType>(
        "Oil feed", "Where oil enters the film.", "geom.feed_type", case_.feed.type,
        {{gui::FeedType::Hole, "Hole", "Circular feed hole in the bearing."},
         {gui::FeedType::Groove, "Groove", "Axial feed groove across the bearing length."}});
    changed |= ui::form_number("Feed angle",
                               "Centre of the feed, measured from the fixed-frame x axis in the direction "
                               "of rotation (theta_in = 0 at the feed).",
                               "geom.feed_angle", case_.feed.angle_deg, gui::UnitFamily::Angle, units());
    if (case_.feed.type == gui::FeedType::Hole) {
        changed |= ui::form_number("Hole diameter", "Diameter of the feed hole. Typical: 1 to 4 mm.",
                                   "geom.hole_diameter", case_.feed.hole_diameter, gui::UnitFamily::Length,
                                   units());
    } else {
        changed |= ui::form_number("Groove width", "Angular width of the axial feed groove. Typical: 10 to 30 deg.",
                                   "geom.groove_width", case_.feed.groove_width_deg, gui::UnitFamily::Angle,
                                   units());
    }

    if (advanced) {
        if (case_.feed.type == gui::FeedType::Hole) {
            bool mid_length = case_.feed.hole_axial_position < 0.0;
            if (ui::form_checkbox("Hole at mid-length", "Feed hole centred along the bearing length.",
                                  "geom.hole_mid", mid_length)) {
                case_.feed.hole_axial_position = mid_length ? -1.0 : 0.5 * case_.bearing_length;
                changed = true;
            }
            if (!mid_length) {
                changed |= ui::form_number("Hole axial position",
                                           "Distance of the hole centre from the bearing end at z = 0.",
                                           "geom.hole_axial", case_.feed.hole_axial_position,
                                           gui::UnitFamily::Length, units());
            }
        }
        changed |= ui::form_number("Tilt angle gamma",
                                   "Shaft axis relative to the bearing axis (slope of the film along the "
                                   "length). 0 = aligned. Typical: below 1 mrad.",
                                   "geom.tilt", case_.tilt_angle, gui::UnitFamily::TiltAngle, units());
        changed |= ui::form_number("Tilt direction",
                                   "Fixed-frame angle toward which the shaft axis tilts, in the direction "
                                   "of rotation.",
                                   "geom.tilt_direction", case_.tilt_direction_deg, gui::UnitFamily::Angle,
                                   units());
    }

    char text[64];
    if (case_.bearing_diameter > 0.0) {
        std::snprintf(text, sizeof(text), "%.3g", case_.bearing_length / case_.bearing_diameter);
        ui::form_value("L/D", "Length-to-diameter ratio.", text);
        std::snprintf(text, sizeof(text), "%.3g", case_.radial_clearance / (0.5 * case_.bearing_diameter));
        ui::form_value("c/R", "Relative clearance: radial clearance per radius.", text);
    }
    ui::end_form();
    if (changed) mark_changed();
}
