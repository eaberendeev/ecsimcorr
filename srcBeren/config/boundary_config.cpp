#include "boundary_config.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "json_settings.h"
#include "runtime/profiler.h"
#include "simulation/boundary_conditions.h"

namespace {
static inline Face string_to_face(const std::string& s) {
    std::string upper = s;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return std::toupper(c); });

    if (upper == "XMIN")
        return Face::XMIN;
    if (upper == "XMAX")
        return Face::XMAX;
    if (upper == "YMIN")
        return Face::YMIN;
    if (upper == "YMAX")
        return Face::YMAX;
    if (upper == "ZMIN")
        return Face::ZMIN;
    if (upper == "ZMAX")
        return Face::ZMAX;
    if (upper == "CYLINDER")
        return Face::CYLINDER;

    // Если ничего не подошло — бросаем исключение
    throw std::invalid_argument("Unknown face string: " + s);
}

SecondaryEmissionModel secondary_emission_from_config(Face face, const nlohmann::json& params) {
    const std::string product = params.value("product", "Electrons");
    if (!params.contains("sources") || !params["sources"].is_array() || params["sources"].empty()) {
        throw std::runtime_error(
            R"msg("second_emission" requires a non-empty "sources" array, e.g.:
  {"second_emission": {"face": "ZMIN", "product": "Electrons", "sources": [
      {"species": "Ions", "yield": 0.3, "energy": {"type": "fixed", "kev": 2.0}},
      {"species": "Electrons", "yield": 0.1, "energy": {"type": "temperature", "temperature_kev": [1.0,1.0,1.0]}},
      {"species": "Ions2", "yield": 0.5, "energy": {"type": "fraction", "fraction": 0.5}}
    ]}})msg");
    }
    std::vector<EmissionSourceRule> rules;
    for (const auto& src : params["sources"]) {
        EmissionSourceRule rule;
        rule.species = src.at("species").get<std::string>();
        if (src.contains("yield") && src["yield"].is_number()) {
            rule.yield_model = EmissionSourceRule::YieldModel::Constant;
            rule.yield = src["yield"].get<double>();
            if (rule.yield < 0.0) {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": \"yield\" must be >= 0, got " + std::to_string(rule.yield));
            }
        } else if (src.contains("yield") && src["yield"].is_object()) {
            const nlohmann::json& ycfg = src["yield"];
            const std::string ytype = ycfg.value("type", "");
            if (ytype == "vaughan") {
                rule.yield_model = EmissionSourceRule::YieldModel::Vaughan;
                rule.delta_max = ycfg.at("delta_max").get<double>();
                rule.energy_max_kev = ycfg.at("energy_max_kev").get<double>();
                rule.threshold_kev = ycfg.value("threshold_kev", 0.0);
                if (rule.delta_max <= 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"delta_max\" must be > 0, got " + std::to_string(rule.delta_max));
                }
                if (rule.energy_max_kev <= 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"energy_max_kev\" must be > 0, got " +
                                             std::to_string(rule.energy_max_kev));
                }
                if (rule.threshold_kev < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"threshold_kev\" must be >= 0, got " +
                                             std::to_string(rule.threshold_kev));
                }
            } else if (ytype == "threshold") {
                rule.yield_model = EmissionSourceRule::YieldModel::Threshold;
                rule.yield = ycfg.at("yield").get<double>();
                rule.threshold_kev = ycfg.at("threshold_kev").get<double>();
                if (rule.yield < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"yield\" must be >= 0, got " + std::to_string(rule.yield));
                }
                if (rule.threshold_kev < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"threshold_kev\" must be >= 0, got " +
                                             std::to_string(rule.threshold_kev));
                }
            } else {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": unknown \"yield\" type \"" + ytype +
                                         "\" (expected number, {\"type\":\"vaughan\",...} or "
                                         "{\"type\":\"threshold\",...})");
            }
        } else {
            throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                     "\": \"yield\" must be a number or an object");
        }
        if (!src.contains("energy") || !src["energy"].is_object()) {
            throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                     "\" requires an \"energy\" object, e.g. {\"type\": \"fixed\", \"kev\": 2.0}");
        }
        const auto& energy_cfg = src["energy"];
        const std::string etype = energy_cfg.at("type").get<std::string>();
        if (etype == "fixed") {
            rule.energy_type = EmissionSourceRule::EnergyType::Fixed;
            rule.fixed_kev = energy_cfg.at("kev").get<double>();
            if (rule.fixed_kev < 0.0) {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": \"kev\" must be >= 0, got " + std::to_string(rule.fixed_kev));
            }
        } else if (etype == "temperature") {
            rule.energy_type = EmissionSourceRule::EnergyType::Temperature;
            rule.temperature_kev = util::parse_double3(energy_cfg.at("temperature_kev"));
            if (rule.temperature_kev.x() < 0.0 || rule.temperature_kev.y() < 0.0 || rule.temperature_kev.z() < 0.0) {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": every component of \"temperature_kev\" must be >= 0, got [" +
                                         std::to_string(rule.temperature_kev.x()) + ", " +
                                         std::to_string(rule.temperature_kev.y()) + ", " +
                                         std::to_string(rule.temperature_kev.z()) + "]");
            }
            rule.temperature_mean =
                energy_cfg.contains("mean") ? util::parse_double3(energy_cfg["mean"]) : Vector3R(0, 0, 0);
        } else if (etype == "fraction") {
            rule.energy_type = EmissionSourceRule::EnergyType::Fraction;
            rule.fraction = energy_cfg.at("fraction").get<double>();
            if (rule.fraction < 0.0 || rule.fraction > 1.0) {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": \"fraction\" must be in [0..1], got " + std::to_string(rule.fraction));
            }
        } else {
            throw std::runtime_error("\"second_emission\" source \"" + rule.species + "\": unknown energy type \"" +
                                     etype + "\" (expected \"fixed\", \"temperature\" or \"fraction\")");
        }
        rules.push_back(rule);
    }
    for (const auto& rule : rules) {
        if (rule.yield_model == EmissionSourceRule::YieldModel::Constant && rule.species == product &&
            rule.yield >= 1.0) {
            throw std::runtime_error("\"second_emission\": source \"" + rule.species + "\" equals product \"" +
                                     product +
                                     "\" with constant yield >= 1 — self-amplifying emission would run away; use "
                                     "\"vaughan\"/\"threshold\" with an energy gate or yield < 1");
        }
    }
    return SecondaryEmissionModel(face, product, rules);
}

void add_boundary_condition(std::vector<std::unique_ptr<FieldBoundaryCondition>>& conditions_,
                            BoundaryEmissionModels& emissions_, std::array<bool, 3>& periodic_, const std::string& type,
                            const nlohmann::json& params, const Domain& domain) {
    const std::string face_str = params.at("face");
    const Face face = string_to_face(face_str);

    if (face == Face::CYLINDER && !domain.geom.use_cylinder) {
        throw std::runtime_error(
            R"msg(Boundary condition uses face "CYLINDER", but "CylinderDomain" is not configured.
Add to system_config:
  "CylinderDomain": {"radius": <value>, "center": [<x>, <y>]})msg");
    }

    if (type == "bphi") {
        const double electron_threshold_energy = params.at("electron_threshold_energy_");
        const double radius = params.at("radius");
        const double gap = params.at("gap");
        conditions_.push_back(
            std::make_unique<ElectronReflectionCondition>(face, radius, electron_threshold_energy, gap));
    } else if (type == "second_emission") {
        // ВАЖНО: эмиссия — НЕ consuming-условие: ничего не добавляем в
        // conditions_, только в emissions_.
        emissions_[static_cast<size_t>(face)].push_back(secondary_emission_from_config(face, params));
    } else if (type == "second_emisson") {
        // Устаревший вариант написания (опечатка): сообщаем и просим переименовать.
        throw std::runtime_error(
            R"msg(Boundary condition "second_emisson" is a typo and no longer supported.
Rename it to "second_emission" and use the new format, e.g.:
  {"second_emission": {"face": "ZMIN", "product": "Electrons", "sources": [
      {"species": "Ions", "yield": 0.3, "energy": {"type": "fixed", "kev": 2.0}}]}})msg");
    } else if (type == "periodic") {
        // Periodic folding of fields/operators acts on the whole axis at once
        // (it pairs the low and high boundary layers mutually), so a single
        // condition per axis is sufficient. XMIN+XMAX entries in the config
        // describe the same axis; adding both would apply the fold twice and
        // double the field values in the boundary layers.
        int axis = -1;
        if (face == Face::XMIN || face == Face::XMAX) {
            axis = 0;
        } else if (face == Face::YMIN || face == Face::YMAX) {
            axis = 1;
        } else if (face == Face::ZMIN || face == Face::ZMAX) {
            axis = 2;
        }
        if (axis >= 0 && periodic_[axis]) {
            std::cout << "BoundaryConditionHandler: duplicate periodic condition on face " << face_str
                      << " is ignored (axis already periodic)\n";
            return;
        }
        conditions_.push_back(std::make_unique<PeriodicBoundaryCondition>(face));
        if (axis >= 0) {
            periodic_[axis] = true;
        }
    } else if (type == "electron_reflection") {
        const double radius = params.value("radius", 5.0);
        const double energy_threshold = params.value("energy_threshold", 0.1) / SGS::MC2;
        conditions_.push_back(std::make_unique<ElectronReflectionCondition>(face, radius, energy_threshold));
    } else if (type == "er0") {
        if (face != Face::ZMIN && face != Face::ZMAX) {
            throw std::runtime_error("\"er0\": \"face\" must be \"ZMIN\" or \"ZMAX\", got \"" + face_str + "\"");
        }
        const std::string profile_name = params.value("profile", "step");
        Er0Profile profile = Er0Profile::Step;
        if (profile_name == "step") {
            profile = Er0Profile::Step;
        } else if (profile_name == "cosine" || profile_name == "cos") {
            profile = Er0Profile::Cosine;
        } else {
            throw std::runtime_error("\"er0\": unsupported \"profile\": \"" + profile_name +
                                     "\" (supported values: \"step\", \"cosine\")");
        }
        const double inner_radius = params.value("inner_radius", 5.0);
        // Для косинусного профиля r0 — центр кольца; по умолчанию совпадает с inner_radius.
        const double r0 = params.value("r0", inner_radius);
        const double width = params.value("width", 2.0);
        const double potential_drop = params.value("potential_drop", 0.1) / SGS::MC2;
        if (!std::isfinite(inner_radius) || !(inner_radius >= 0.0)) {
            throw std::runtime_error("\"er0\": \"inner_radius\" must be >= 0, got " + std::to_string(inner_radius));
        }
        if (!std::isfinite(r0) || !(r0 >= 0.0)) {
            throw std::runtime_error("\"er0\": \"r0\" must be >= 0, got " + std::to_string(r0));
        }
        if (!std::isfinite(width) || !(width > 0.0)) {
            throw std::runtime_error("\"er0\": \"width\" must be > 0, got " + std::to_string(width));
        }
        if (!std::isfinite(potential_drop)) {
            throw std::runtime_error("\"er0\": \"potential_drop\" must be finite, got " +
                                     std::to_string(potential_drop));
        }
        if (profile == Er0Profile::Cosine && r0 < 0.5 * width) {
            throw std::runtime_error("\"er0\" cosine profile requires r0 >= width/2 for a full potential drop");
        }
        const std::string mode = params.value("mode", "overwrite");
        if (mode != "overwrite" && mode != "solv") {
            throw std::runtime_error("\"er0\": unsupported \"mode\": \"" + mode +
                                     "\" (supported values: \"overwrite\", \"solv\")");
        }
        conditions_.push_back(
            std::make_unique<Er0Condition>(face, profile, inner_radius, r0, width, potential_drop, mode));
    } else {
        std::cerr << "Unknown boundary condition type: " << type << std::endl;
    }
}

}   // namespace

BoundaryConditionHandler boundary_conditions_from_config(const nlohmann::json& sys_config, const Domain& domain) {
    RECORD_TIMER;
    std::vector<std::unique_ptr<FieldBoundaryCondition>> conditions_;
    BoundaryEmissionModels emissions_;
    std::array<bool, 3> periodic_{};
    if (!sys_config.contains("Boundary_conditions")) {
        return BoundaryConditionHandler(std::move(conditions_), std::move(emissions_), periodic_);
    }
    const nlohmann::json& config = sys_config["Boundary_conditions"];
    if (!config.is_array()) {
        std::cerr << "BoundaryConditionHandler: expected array in JSON\n";
        return BoundaryConditionHandler(std::move(conditions_), std::move(emissions_), periodic_);
    }

    for (const nlohmann::json& item : config) {
        if (item.empty())
            continue;

        nlohmann::json::const_iterator it = item.begin();
        const std::string type = it.key();
        const nlohmann::json& params = it.value();

        // "open" accepts only an array of faces; all faces are gathered into
        // a single OpenBoundaryConditionArray (this class exists exactly for
        // such batching)
        if (type == "open") {
            if (!params.is_array()) {
                throw std::runtime_error(
                    R"msg("open" expects an array of faces, e.g. "open": [{"face": "XMIN"}, {"face": "XMAX"}])msg");
            }
            std::vector<Face> faces;
            for (const nlohmann::json& faceJs : params) {
                faces.push_back(string_to_face(faceJs.at("face")));
            }
            conditions_.push_back(std::make_unique<OpenBoundaryConditionArray>(faces));
            continue;
        }

        // Other types accept either a single face object or an array of them
        // (e.g. {"periodic": [{"face": "YMIN"}, {"face": "YMAX"}]})
        if (params.is_array()) {
            for (const nlohmann::json& faceJs : params) {
                add_boundary_condition(conditions_, emissions_, periodic_, type, faceJs, domain);
            }
        } else {
            add_boundary_condition(conditions_, emissions_, periodic_, type, params, domain);
        }
    }

    // "er0" заменяет "open" на своей грани: дублировать open нельзя — в solv
    // это двусмысленно (два разных маскирования), в overwrite избыточно.
    for (const auto& cond : conditions_) {
        const auto* er0 = dynamic_cast<const Er0Condition*>(cond.get());
        if (!er0)
            continue;
        const Face er0_face = er0->face();
        for (const auto& other : conditions_) {
            if (other.get() == cond.get())
                continue;
            if (dynamic_cast<const Er0Condition*>(other.get()))
                continue;
            bool covers = false;
            if (const auto* open = dynamic_cast<const OpenBoundaryCondition*>(other.get())) {
                covers = open->face() == er0_face;
            } else if (const auto* arr = dynamic_cast<const OpenBoundaryConditionArray*>(other.get())) {
                for (int i = 0; i < arr->createdBc_ && !covers; ++i) covers = arr->faces_[i] == er0_face;
            }
            if (covers) {
                throw std::runtime_error("\"er0\" on face " + std::to_string(static_cast<int>(er0_face)) +
                                         " replaces \"open\": remove the duplicate \"open\" entry for the same face");
            }
        }
    }

    // Проверка ПОСЛЕ обработки всех записей: для каждой грани, на которой
    // задана вторичная эмиссия, должно существовать consuming-условие (open,
    // bphi и т.п.) на той же грани — иначе эмиссия никогда не сработает
    // (частицы не будут удалены на этой грани). Условие может идти в JSON
    // после second_emission, поэтому проверяем только в конце цикла.
    for (size_t f = 0; f < emissions_.size(); ++f) {
        const Face face = static_cast<Face>(f);
        const auto& models = emissions_[f];
        if (models.empty())
            continue;
        bool has_consumer = false;
        for (const auto& cond : conditions_) {
            if (cond->face() == face) {
                has_consumer = true;
                break;
            }
            // OpenBoundaryConditionArray покрывает сразу несколько граней.
            if (const auto* arr = dynamic_cast<const OpenBoundaryConditionArray*>(cond.get())) {
                for (int i = 0; i < arr->createdBc_; ++i) {
                    if (arr->faces_[i] == face) {
                        has_consumer = true;
                        break;
                    }
                }
            }
            if (has_consumer)
                break;
        }
        if (!has_consumer) {
            throw std::runtime_error("\"second_emission\" on face " + std::to_string(static_cast<int>(face)) +
                                     " needs a consuming boundary condition (e.g. \"open\" or \"bphi\") on the same "
                                     "face, otherwise it never fires.");
        }
    }

    return BoundaryConditionHandler(std::move(conditions_), std::move(emissions_), periodic_);
}
