#ifndef CECE_DIAGNOSTICS_HPP
#define CECE_DIAGNOSTICS_HPP

#include <string>
#include <unordered_map>
#include <vector>

#include "cece/cece_compute.hpp"
#include "cece/cece_config.hpp"

namespace cece {

/**
 * @struct DiagnosticDeclaration
 * @brief A diagnostic recorded before its field can be allocated.
 *
 * Schemes take this path during Phase 1 (no grid known yet, or no manager
 * bound). Only metadata is retained; the requested placeholder dimensions are
 * intentionally ignored once the manager grid is set (FR-008).
 */
struct DiagnosticDeclaration {
    std::string name;       ///< Unique diagnostic name.
    std::string units;      ///< CF-compliant units (may be empty).
    std::string long_name;  ///< Descriptive name (may be empty).
};

/**
 * @struct DiagnosticFieldInfo
 * @brief A registered diagnostic: a grid-sized field plus its metadata.
 *
 * The DualView3D is shared process-wide for a given name (FR-003): every
 * scheme that registers the same name receives a handle aliasing this storage.
 */
struct DiagnosticFieldInfo {
    DualView3D data;        ///< Grid-sized field (extents == manager grid).
    std::string units;      ///< CF-compliant units.
    std::string long_name;  ///< Descriptive name.
};

/**
 * @class CeceDiagnosticManager
 * @brief Registry of diagnostic fields, sized to the realized grid.
 *
 * Lifecycle (data-model E3): diagnostics declared before the grid is known are
 * held as pending declarations; SetGrid promotes them (allocating at the real
 * grid extent) and makes all later registrations allocate immediately. The
 * registry is framework-free — no AMIO/ESMF assumptions — so its semantics are
 * unit-testable without a coupled framework (research R8 tier-1).
 */
class CeceDiagnosticManager {
   public:
    CeceDiagnosticManager() = default;

    /**
     * @brief Sets the authoritative grid dimensions and promotes pending fields.
     *
     * After this call every diagnostic allocates at (nx, ny, nz) regardless of
     * the placeholder dims requested at declaration time (FR-008). May be
     * called only once with a given extent; a conflicting re-set throws.
     *
     * @param nx Grid dimension in x.
     * @param ny Grid dimension in y.
     * @param nz Grid dimension in z.
     */
    void SetGrid(int nx, int ny, int nz);

    /**
     * @brief Declares names that diagnostics must never use (species names).
     * @param names Reserved names (typically the configured species).
     */
    void SetReservedNames(std::vector<std::string> names);

    /**
     * @brief Registers a diagnostic, returning the shared grid-sized field.
     *
     * Duplicate names return the existing field (FR-003). Coordinate names and
     * reserved species names throw, naming the offender (FR-009). Metadata is
     * first-non-empty-wins: a later non-empty units/long_name fills an unset
     * value but never clobbers an established one (FR-004, FR-016).
     *
     * Before SetGrid the declaration is recorded as pending and an empty handle
     * is returned; the field is allocated when the grid is set.
     *
     * @return DualView3D the shared field (grid-sized once the grid is set), or
     *         an empty handle while the declaration is still pending.
     */
    DualView3D RegisterDiagnostic(const std::string& name, int nx, int ny, int nz, const std::string& units = "", const std::string& long_name = "");

    /// @brief Read-only view of allocated fields (name -> field + metadata).
    const std::unordered_map<std::string, DiagnosticFieldInfo>& Fields() const {
        return fields_;
    }

    /// @brief Names of allocated fields in stable first-seen order (cap/manifest).
    const std::vector<std::string>& Order() const {
        return order_;
    }

    /// @brief True once SetGrid has established the grid extent.
    bool HasGrid() const {
        return grid_set_;
    }

    int nx() const {
        return nx_;
    }
    int ny() const {
        return ny_;
    }
    int nz() const {
        return nz_;
    }

   private:
    // Rejects coordinate/bounds/mesh names and reserved species names.
    void ValidateName(const std::string& name) const;
    // Allocates a grid-sized, zero-initialized field and records its name.
    DiagnosticFieldInfo& Allocate(const DiagnosticDeclaration& declaration);

    bool grid_set_ = false;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;

    std::vector<std::string> reserved_names_;
    std::unordered_map<std::string, DiagnosticDeclaration> pending_;
    std::unordered_map<std::string, DiagnosticFieldInfo> fields_;
    std::vector<std::string> order_;  // first-seen order over pending_ then fields_
};

}  // namespace cece

#endif  // CECE_DIAGNOSTICS_HPP
