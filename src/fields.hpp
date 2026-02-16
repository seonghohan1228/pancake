#pragma once
#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include "field.hpp"

class Fields {
private:
    // Owns the memory for all fields
    std::map<std::string, std::unique_ptr<Field>> _storage;

public:
    // 1. Create & Store a new field
    // Usage: fields.add("pressure", mesh);
    Field& add(const std::string& name, const Mesh& mesh, int ghost_layers = 2, GridLocation loc = GridLocation::CENTER) {
        if (_storage.find(name) != _storage.end()) {
            throw std::runtime_error("Field already exists: " + name);
        }
        _storage[name] = std::make_unique<Field>(name, mesh, ghost_layers, loc);
        return *_storage[name];
    }

    // 2. Access existing field (Read/Write)
    // Usage: fields["pressure"] = ...
    Field& operator[](const std::string& name) {
        auto it = _storage.find(name);
        if (it == _storage.end()) {
            throw std::runtime_error("Field not found: " + name);
        }
        return *(it->second);
    }

    // 3. Read-only access
    const Field& operator[](const std::string& name) const {
        auto it = _storage.find(name);
        if (it == _storage.end()) {
            throw std::runtime_error("Field not found: " + name);
        }
        return *(it->second);
    }

    // Check if field exists
    bool exists(const std::string& name) const {
        return _storage.find(name) != _storage.end();
    }

    // Iterators (loop: for(auto& [name, field] : fields) ... )
    auto begin() { return _storage.begin(); }
    auto end() { return _storage.end(); }
};
