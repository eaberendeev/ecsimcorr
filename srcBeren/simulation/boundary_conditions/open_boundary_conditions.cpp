#include "runtime/profiler.h"
#include "simulation/boundary_conditions.h"

namespace {
Vector3R electric_position(const Domain& domain, int stored_component) {
    const std::array<int, 4> index = domain.grid.pos_vind_range(stored_component);
    return domain.storage_node_position(index[0], index[1], index[2], FieldType::ELECTRIC, index[3]);
}

// Zero the rows and the columns of the unknowns marked outside.
void zero_rows_and_columns(Operator& mat, const std::vector<char>& outside) {
    mat.makeCompressed();
    double* values = mat.valuePtr();
    const int* outer = mat.outerIndexPtr();
    const int* inner = mat.innerIndexPtr();
#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < mat.rows(); ++row)
        for (int p = outer[row]; p < outer[row + 1]; ++p)
            if (outside[row] || outside[inner[p]])
                values[p] = 0.0;
}
}   // namespace

void OpenBoundaryConditionArray::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;

    auto set_zero = [this](double& value, int i, int j, int k, int d, const Domain& dom, FieldType ft) {
        const Vector3R pos = dom.storage_node_position(i, j, k, ft, d);
        for (int numBc = 0; numBc < createdBc_; ++numBc) {
            if (dom.geom.is_outside_face(faces_[numBc], pos, -gaps_[numBc] + epss_[numBc])) {
                value = 0.0;
                break;
            }
        }
    };

    if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d) set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::ELECTRIC);
    } else if (field_type == FieldType::MAGNETIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d) set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::MAGNETIC);
    } else if (field_type == FieldType::DENSITY) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k) set_zero(field(i, j, k, 0), i, j, k, 0, domain, FieldType::DENSITY);
    } else {
        std::cerr << "OpenBoundaryCondition: unsupported field type\n";
        assert(false);
    }
}

void OpenBoundaryConditionArray::apply_to_operator(Operator& mat, const Domain& domain,
                                                   const std::vector<int>& storage_index) {
    RECORD_TIMER_PARAMS(sizeof(int) * (mat.rows() + mat.nonZeros()), timer::MeasureUnit::byte);
    std::vector<char> outside(mat.rows());
#pragma omp parallel for schedule(static)
    for (int row = 0; row < mat.rows(); ++row) {
        const Vector3R position = electric_position(domain, storage_index[row]);
        bool beyond = false;
        for (int bc = 0; bc < createdBc_ && !beyond; ++bc)
            beyond = domain.geom.is_outside_face(faces_[bc], position, -gaps_[bc] + epss_[bc]);
        outside[row] = beyond;
    }
    zero_rows_and_columns(mat, outside);
}

ParticleFateResult OpenBoundaryConditionArray::classify_particle(const Particle& p, const ParticlesArray& /*particles*/,
                                                                 const Domain& domain, Face hit_face) const {
    for (int i = 0; i < createdBc_; ++i) {
        if (hit_face != Face::UNAVAILABLE ? hit_face == faces_[i] : domain.geom.is_outside_face(faces_[i], p.coord)) {
            return {ParticleFate::Removed, faces_[i]};
        }
    }
    return {};
}

void OpenBoundaryCondition::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;

    auto set_zero = [gap = gap_, eps = eps_](double& value, int i, int j, int k, int d, Face face, const Domain& dom,
                                             FieldType ft) {
        const Vector3R pos = dom.storage_node_position(i, j, k, ft, d);
        if (dom.geom.is_outside_face(face, pos, -gap + eps))
            value = 0.0;
    };

    if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d)
                        set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::ELECTRIC);
    } else if (field_type == FieldType::MAGNETIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d)
                        set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::MAGNETIC);
    } else if (field_type == FieldType::DENSITY) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    set_zero(field(i, j, k, 0), i, j, k, 0, face_, domain, FieldType::DENSITY);
    } else {
        std::cerr << "OpenBoundaryCondition: unsupported field type\n";
        assert(false);
    }
}

void OpenBoundaryCondition::apply_to_operator(Operator& mat, const Domain& domain,
                                              const std::vector<int>& storage_index) {
    RECORD_TIMER_PARAMS(sizeof(int) * (mat.rows() + mat.nonZeros()), timer::MeasureUnit::byte);
    std::vector<char> outside(mat.rows());
#pragma omp parallel for schedule(static)
    for (int row = 0; row < mat.rows(); ++row)
        outside[row] = domain.geom.is_outside_face(face_, electric_position(domain, storage_index[row]), -gap_ + eps_);
    zero_rows_and_columns(mat, outside);
}

ParticleFateResult OpenBoundaryCondition::classify_particle(const Particle& p, const ParticlesArray& /*particles*/,
                                                            const Domain& domain, Face hit_face) const {
    if (hit_face != Face::UNAVAILABLE ? hit_face != face_ : !domain.geom.is_outside_face(face_, p.coord))
        return {};
    return {ParticleFate::Removed, face_};
}
