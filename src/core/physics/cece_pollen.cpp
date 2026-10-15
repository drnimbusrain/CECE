#include "cece/physics/cece_pollen.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "cece/cece_physics_factory.hpp"

namespace cece {
namespace {

KOKKOS_INLINE_FUNCTION
double linear_suppression(double value, double low, double high) {
    if (value <= low) return 1.0;
    if (value >= high) return 0.0;
    return (high - value) / (high - low);
}

}  // namespace

static PhysicsRegistration<PollenScheme> register_scheme("pollen");
static PhysicsRegistration<PollenScheme> register_artemisia_scheme("pollen_artemisia");
static PhysicsRegistration<PollenScheme> register_mugwort_scheme("pollen_mugwort");
static PhysicsRegistration<PollenScheme> register_chenopod_scheme("pollen_chenopod");
static PhysicsRegistration<PollenScheme> register_ragweed_scheme("pollen_ragweed");
static PhysicsRegistration<PollenScheme> register_grass_scheme("pollen_grass");
static PhysicsRegistration<PollenScheme> register_alder_scheme("pollen_alder");
static PhysicsRegistration<PollenScheme> register_ash_scheme("pollen_ash");
static PhysicsRegistration<PollenScheme> register_birch_scheme("pollen_birch");
static PhysicsRegistration<PollenScheme> register_cottonwood_scheme("pollen_cottonwood");
static PhysicsRegistration<PollenScheme> register_cypress_scheme("pollen_cypress");
static PhysicsRegistration<PollenScheme> register_elm_scheme("pollen_elm");
static PhysicsRegistration<PollenScheme> register_hazel_scheme("pollen_hazel");
static PhysicsRegistration<PollenScheme> register_juniper_scheme("pollen_juniper");
static PhysicsRegistration<PollenScheme> register_maple_scheme("pollen_maple");
static PhysicsRegistration<PollenScheme> register_oak_scheme("pollen_oak");
static PhysicsRegistration<PollenScheme> register_olive_scheme("pollen_olive");
static PhysicsRegistration<PollenScheme> register_pine_scheme("pollen_pine");
static PhysicsRegistration<PollenScheme> register_plane_scheme("pollen_plane");
static PhysicsRegistration<PollenScheme> register_nettle_scheme("pollen_nettle");
static PhysicsRegistration<PollenScheme> register_total_scheme("pollen_total");

void PollenScheme::Initialize(const conf::Value& config, CeceDiagnosticManager* diag_manager) {
    BasePhysicsScheme::Initialize(config, diag_manager);

    for (const auto& key : config.keys()) {
        if (key == "phenology_method" || key == "season_start_doy" || key == "season_end_doy" ||
            key == "vegetation_fraction_default" || key == "gaussian_width" || key == "use_autumn_trigger" ||
            key.starts_with("autumn_") || key == "precipitation_low_mm_interval" || key == "precipitation_high_mm_interval") {
            throw std::invalid_argument("PollenScheme: legacy option '" + key + "' is unsupported; migrate to inline base_pollen_potential (docs/pollen.md)");
        }
    }
    for (const auto& key : config["input_mapping"].keys()) {
        if (key == "annual_pollen_production" || key == "vegetation_fraction" || key == "convective_velocity" ||
            key == "precipitation" || key == "day_of_year" || key == "sunshine_hours" ||
            key == "season_start_doy" || key == "season_end_doy" || key == "autumn_accumulated_forcing") {
            throw std::invalid_argument("PollenScheme: legacy input mapping '" + key + "'; see docs/pollen.md");
        }
    }
    if (config["output_mapping"]["pollen_phenology_forcing"]) {
        throw std::invalid_argument("PollenScheme: pollen_phenology_forcing is no longer an output");
    }
    const auto option = [&config](const std::string& key, double fallback) {
        return config[key].is_defined() ? config[key].as_double() : fallback;
    };
    roughness_length_m_ = option("roughness_length_m", roughness_length_m_);
    friction_velocity_threshold_m_s_ = option("friction_velocity_threshold_m_s", friction_velocity_threshold_m_s_);
    friction_velocity_scale_m_s_ = option("friction_velocity_scale_m_s", friction_velocity_scale_m_s_);
    precipitation_low_mm_h_ = option("precipitation_low_mm_h", precipitation_low_mm_h_);
    precipitation_high_mm_h_ = option("precipitation_high_mm_h", precipitation_high_mm_h_);
    rh_low_percent_ = option("rh_low_percent", rh_low_percent_);
    rh_high_percent_ = option("rh_high_percent", rh_high_percent_);
    temperature_threshold_c_ = option("temperature_threshold_c", temperature_threshold_c_);
    temperature_slope_ = option("temperature_slope", temperature_slope_);
    particle_diameter_m_ = option("particle_diameter_um", 30.0) * 1.0e-6;
    particle_density_kg_m3_ = option("particle_density_kg_m3", particle_density_kg_m3_);
    for (double value : {roughness_length_m_, friction_velocity_threshold_m_s_, friction_velocity_scale_m_s_,
                         precipitation_low_mm_h_, precipitation_high_mm_h_, rh_low_percent_, rh_high_percent_,
                         temperature_threshold_c_, temperature_slope_, particle_diameter_m_, particle_density_kg_m3_}) {
        if (!std::isfinite(value)) throw std::invalid_argument("PollenScheme: options must be finite");
    }
    if (roughness_length_m_ <= 0.0 || roughness_length_m_ >= 10.0 ||
        friction_velocity_threshold_m_s_ < 0.0 || friction_velocity_scale_m_s_ <= 0.0) {
        throw std::invalid_argument("PollenScheme: require 0 < roughness_length_m < 10, nonnegative friction threshold and positive friction scale");
    }
    if (precipitation_low_mm_h_ < 0.0 || precipitation_high_mm_h_ <= precipitation_low_mm_h_ ||
        rh_low_percent_ < 0.0 || rh_high_percent_ > 100.0 || rh_high_percent_ <= rh_low_percent_) {
        throw std::invalid_argument("PollenScheme: require increasing nonnegative rain thresholds and 0 <= RH low < RH high <= 100");
    }
    if (particle_diameter_m_ <= 0.0 || particle_density_kg_m3_ <= 0.0 || temperature_slope_ <= 0.0) {
        throw std::invalid_argument("PollenScheme: particle properties and temperature_slope must be positive");
    }
}

void PollenScheme::Run(CeceImportState& import_state, CeceExportState& export_state) {
    auto potential = ResolveImport("base_pollen_potential", import_state);
    auto wind_speed = ResolveImport("wind_speed", import_state);
    auto precipitation = ResolveImport("precipitation_rate", import_state);
    auto relative_humidity = ResolveImport("relative_humidity", import_state);
    auto temperature = ResolveImport("temperature", import_state);

    auto number_flux = ResolveExport("pollen_number_emissions", export_state);
    auto mass_flux = ResolveExport("pollen_mass_emissions", export_state);
    auto diameter = ResolveExport("pollen_diameter", export_state);
    auto density = ResolveExport("pollen_density", export_state);
    RequireFields("PollenScheme::Run", {{MapInput("base_pollen_potential"), potential.data() != nullptr},
                                        {MapInput("wind_speed"), wind_speed.data() != nullptr},
                                        {MapInput("precipitation_rate"), precipitation.data() != nullptr},
                                        {MapInput("relative_humidity"), relative_humidity.data() != nullptr},
                                        {MapInput("temperature"), temperature.data() != nullptr},
                                        {MapOutput("pollen_number_emissions"), number_flux.data() != nullptr},
                                        {MapOutput("pollen_mass_emissions"), mass_flux.data() != nullptr}});
    const int nx = static_cast<int>(number_flux.extent(0));
    const int ny = static_cast<int>(number_flux.extent(1));
    const auto check_shape = [nx, ny](const auto& field) {
        if (field.extent(0) != static_cast<std::size_t>(nx) || field.extent(1) != static_cast<std::size_t>(ny) || field.extent(2) < 1) {
            throw std::runtime_error("PollenScheme::Run: fields must share horizontal extents and have a surface layer");
        }
    };
    if (nx <= 0 || ny <= 0) throw std::runtime_error("PollenScheme::Run: empty horizontal domain");
    for (const auto& field : {potential, wind_speed, precipitation, relative_humidity, temperature}) check_shape(field);
    check_shape(number_flux);
    check_shape(mass_flux);
    if (diameter.data()) check_shape(diameter);
    if (density.data()) check_shape(density);

    int invalid_cells = 0;
    Kokkos::parallel_reduce(
        "Pollen_Validate", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
        KOKKOS_LAMBDA(int i, int j, int& invalid) {
            const double base = potential(i, j, 0);
            const double wind = wind_speed(i, j, 0);
            const double rain = precipitation(i, j, 0);
            const double rh = relative_humidity(i, j, 0);
            const double t = temperature(i, j, 0);
            if (!Kokkos::isfinite(base) || base < 0.0 || !Kokkos::isfinite(wind) || wind < 0.0 ||
                !Kokkos::isfinite(rain) || rain < 0.0 || !Kokkos::isfinite(rh) || rh < 0.0 || rh > 100.0 ||
                !Kokkos::isfinite(t) || t <= 0.0) ++invalid;
        }, invalid_cells);
    if (invalid_cells) {
        throw std::runtime_error("PollenScheme::Run: nonfinite or out-of-range surface inputs; require base >= 0, U10 >= 0, rain >= 0, RH in [0,100], T > 0 K");
    }

    const double wind_conversion = 0.4 / std::log(10.0 / roughness_length_m_);
    const double wind_threshold = friction_velocity_threshold_m_s_;
    const double wind_scale = friction_velocity_scale_m_s_;
    const double precipitation_low = precipitation_low_mm_h_;
    const double precipitation_high = precipitation_high_mm_h_;
    const double rh_low = rh_low_percent_;
    const double rh_high = rh_high_percent_;
    const double temperature_threshold = temperature_threshold_c_;
    const double temperature_slope = temperature_slope_;
    const double diameter_m = particle_diameter_m_;
    const double density_kg_m3 = particle_density_kg_m3_;
    const double grain_mass = std::numbers::pi / 6.0 * diameter_m * diameter_m * diameter_m * density_kg_m3;
    const bool has_diameter = diameter.data() != nullptr;
    const bool has_density = density.data() != nullptr;

    Kokkos::parallel_for(
        "Pollen_Emissions", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
            const double temperature_c = temperature(i, j, 0) - 273.15;
            const double friction_velocity = wind_conversion * wind_speed(i, j, 0);
            const double wind_factor = friction_velocity <= wind_threshold ? 0.0 :
                1.0 - std::exp(-(friction_velocity - wind_threshold) / wind_scale);
            const double rain_factor = linear_suppression(precipitation(i, j, 0), precipitation_low, precipitation_high);
            const double humidity_factor = linear_suppression(relative_humidity(i, j, 0), rh_low, rh_high);
            const double temperature_factor = 1.0 / (1.0 + std::exp(-temperature_slope * (temperature_c - temperature_threshold)));
            const double flux = potential(i, j, 0) * wind_factor * rain_factor * humidity_factor * temperature_factor;

            number_flux(i, j, 0) = flux;
            mass_flux(i, j, 0) = flux * grain_mass;
            if (has_diameter) diameter(i, j, 0) = diameter_m;
            if (has_density) density(i, j, 0) = density_kg_m3;
        });

    MarkModified("pollen_number_emissions", export_state);
    MarkModified("pollen_mass_emissions", export_state);
    if (has_diameter) MarkModified("pollen_diameter", export_state);
    if (has_density) MarkModified("pollen_density", export_state);
    Kokkos::fence();
}

}  // namespace cece
