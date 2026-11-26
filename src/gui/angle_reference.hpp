#pragma once

// Selectable zero of every angle axis and polar plot (deck p18):
//   fixed frame, maximum film thickness (theta), minimum film thickness,
//   feed (theta_in), opposite the feed.
// Angles grow in the direction of rotation. Moving zeros (h_max, h_min) are
// taken from the field shown, so a transient view follows the instant.
// Header-only, no toolkit dependency: unit-tested.

#include <cmath>
#include <string>

namespace gui {

enum class AngleZero { Fixed, MaxFilm, MinFilm, Feed, OppositeFeed };

inline const char* angle_zero_name(AngleZero zero) {
    switch (zero) {
        case AngleZero::Fixed: return "Fixed frame";
        case AngleZero::MaxFilm: return "Maximum film (theta)";
        case AngleZero::MinFilm: return "Minimum film";
        case AngleZero::Feed: return "Feed (theta_in)";
        case AngleZero::OppositeFeed: return "Opposite the feed";
    }
    return "Fixed frame";
}

/// Axis label naming the zero, e.g. "Angle from feed, theta_in [deg]".
inline std::string angle_axis_label(AngleZero zero) {
    switch (zero) {
        case AngleZero::Fixed: return "Angle, fixed frame [deg]";
        case AngleZero::MaxFilm: return "Angle from max. film, theta [deg]";
        case AngleZero::MinFilm: return "Angle from min. film [deg]";
        case AngleZero::Feed: return "Angle from feed, theta_in [deg]";
        case AngleZero::OppositeFeed: return "Angle from opposite the feed [deg]";
    }
    return "Angle [deg]";
}

/// Default first value of the axis window (deck p47 shows theta_in from -240).
inline double angle_axis_start(AngleZero zero) {
    switch (zero) {
        case AngleZero::Feed: return -240.0;
        case AngleZero::MinFilm:
        case AngleZero::OppositeFeed: return -180.0;
        default: return 0.0;
    }
}

/// Fixed-frame angle [deg] of the chosen zero.
inline double angle_zero_position(AngleZero zero, double feed_deg, double theta_hmax_deg,
                                  double theta_hmin_deg) {
    switch (zero) {
        case AngleZero::Fixed: return 0.0;
        case AngleZero::MaxFilm: return theta_hmax_deg;
        case AngleZero::MinFilm: return theta_hmin_deg;
        case AngleZero::Feed: return feed_deg;
        case AngleZero::OppositeFeed: return feed_deg + 180.0;
    }
    return 0.0;
}

/// Maps a fixed-frame angle to the displayed angle in [start, start + 360).
/// `rotation_sign` is +1 when the shaft turns toward increasing fixed-frame angle.
inline double display_angle(double fixed_deg, double zero_deg, double start_deg,
                            double rotation_sign = 1.0) {
    double value = (rotation_sign >= 0.0 ? 1.0 : -1.0) * (fixed_deg - zero_deg);
    value = std::fmod(value - start_deg, 360.0);
    if (value < 0.0) value += 360.0;
    return value + start_deg;
}

}  // namespace gui
