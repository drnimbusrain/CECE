#ifndef CECE_POLLEN_HPP
#define CECE_POLLEN_HPP

#include "cece/physics_scheme.hpp"

namespace cece {

class PollenScheme : public BasePhysicsScheme {
   public:
    PollenScheme() = default;
    ~PollenScheme() override = default;

    void Initialize(const conf::Value& config, CeceDiagnosticManager* diag_manager) override;
    void Run(CeceImportState& import_state, CeceExportState& export_state) override;

   private:
    double roughness_length_m_ = 0.1;
    double friction_velocity_threshold_m_s_ = 0.2;
    double friction_velocity_scale_m_s_ = 0.3;
    double precipitation_low_mm_h_ = 0.0;
    double precipitation_high_mm_h_ = 1.0;
    double rh_low_percent_ = 50.0;
    double rh_high_percent_ = 90.0;
    double temperature_threshold_c_ = 10.0;
    double temperature_slope_ = 0.35;
    double particle_diameter_m_ = 30.0e-6;
    double particle_density_kg_m3_ = 1000.0;
};

}  // namespace cece

#endif  // CECE_POLLEN_HPP
