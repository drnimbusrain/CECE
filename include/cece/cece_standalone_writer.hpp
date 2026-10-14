#ifndef CECE_STANDALONE_WRITER_HPP
#define CECE_STANDALONE_WRITER_HPP

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "cece_compute.hpp"
#include "cece_config.hpp"

namespace cece {

/**
 * @struct DiagnosticOut
 * @brief One diagnostic handed to the writer for a single time step.
 *
 * The host view is already synchronized by the caller (the registry owner),
 * so the writer never touches device memory. Metadata comes from the
 * diagnostic registry, which is the single source of truth for units and
 * long_name (contracts/diagnostics-c-api.md §3).
 */
struct DiagnosticOut {
    Kokkos::View<const double***, Kokkos::LayoutLeft, Kokkos::HostSpace> host;  ///< Host-side field values.
    std::string units;                                                          ///< CF units (may be empty -> omitted).
    std::string long_name;                                                      ///< Descriptive name (may be empty -> omitted).
};

class CeceStandaloneWriter {
   public:
    explicit CeceStandaloneWriter(const CeceOutputConfig& config, MPI_Comm comm = MPI_COMM_SELF);
    ~CeceStandaloneWriter();

    CeceStandaloneWriter(const CeceStandaloneWriter&) = delete;
    CeceStandaloneWriter& operator=(const CeceStandaloneWriter&) = delete;
    CeceStandaloneWriter(CeceStandaloneWriter&&) = delete;
    CeceStandaloneWriter& operator=(CeceStandaloneWriter&&) = delete;

    int Initialize(const std::string& start_time_iso8601, int nx, int ny, int nz);

    int InitializeWithCoords(const std::string& start_time_iso8601, int nx, int ny, int nz, const std::vector<double>& lon_coords,
                             const std::vector<double>& lat_coords, const std::string& gridspec_file = "");

    int WriteTimeStep(const std::unordered_map<std::string, DualView3D>& export_fields, const std::map<std::string, DiagnosticOut>& diagnostics,
                      double time_seconds_since_start, int step_index);

    void Finalize();

    bool IsInitialized() const {
        return initialized_;
    }

    /**
     * @brief Seconds since the last diagnostic publication, for cadence.
     * @details The diagnostic publish interval is independent of the file
     * write interval: diagnostics ride file writes (R5) but only when due.
     * -1 means "never published", so the first file write is always due.
     */
    double LastDiagnosticPublishTime() const {
        return last_diag_publish_time_;
    }

    /// @brief Records a diagnostic publication at the given elapsed time.
    void MarkDiagnosticsPublished(double time_seconds_since_start) {
        last_diag_publish_time_ = time_seconds_since_start;
    }

   private:
    CeceOutputConfig config_;
    bool initialized_ = false;
    int record_count_ = 0;
    int nx_ = 0, ny_ = 0, nz_ = 0;
    std::string start_time_iso8601_;
    std::vector<double> lon_coords_;
    std::vector<double> lat_coords_;
    bool use_custom_coords_ = false;
    MPI_Comm comm_ = MPI_COMM_SELF;
    std::string gridspec_file_;
    double last_diag_publish_time_ = -1.0;

    std::string ResolveFilename(double time_seconds_since_start) const;
};

}  // namespace cece

#endif
