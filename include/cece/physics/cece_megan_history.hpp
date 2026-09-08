#ifndef CECE_MEGAN_HISTORY_HPP
#define CECE_MEGAN_HISTORY_HPP

#include <cmath>
#include <conf/value.hpp>
#include <stdexcept>

namespace cece {

/// Scalar effective histories for native C++ MEGAN and MEGAN3.
/// Defaults preserve the previous behavior; no history evolution is performed.
struct MeganHistory {
    double temperature_k = 297.0;
    double par_wm2 = 400.0;
    double days_between_lai = 30.0;
    int day_of_year = 180;
    bool leaf_age_uses_history = false;

    static MeganHistory FromConfig(const conf::Value& config) {
        MeganHistory result;
        if (config["temperature_history_k"].is_defined()) result.temperature_k = config["temperature_history_k"].as_double();
        if (config["par_history_wm2"].is_defined()) result.par_wm2 = config["par_history_wm2"].as_double();
        if (config["days_between_lai"].is_defined()) result.days_between_lai = config["days_between_lai"].as_double();
        if (config["day_of_year"].is_defined()) result.day_of_year = config["day_of_year"].as_int();
        if (config["leaf_age_uses_temperature_history"].is_defined()) {
            result.leaf_age_uses_history = config["leaf_age_uses_temperature_history"].as_bool();
        }
        if (!std::isfinite(result.temperature_k) || result.temperature_k <= 0.0 || !std::isfinite(result.par_wm2) || result.par_wm2 < 0.0 ||
            !std::isfinite(result.days_between_lai) || result.days_between_lai <= 0.0 || result.day_of_year < 1 || result.day_of_year > 366) {
            throw std::invalid_argument(
                "MEGAN history requires positive finite temperature_history_k and days_between_lai, "
                "nonnegative finite par_history_wm2, and day_of_year in [1, 366]");
        }
        return result;
    }
};

}  // namespace cece
#endif
