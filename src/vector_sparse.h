#pragma once

#include <vector>

class VectorSparse {
private:
    struct Element {
        unsigned index;
        float value;
        Element(const unsigned index, const float value) : index(index), value(value) {}
        Element operator-() const {return Element(index, -value);}
    };
    std::vector<Element> data;
public:
    std::size_t size() const { return data.size(); }

    Element& operator[](const std::size_t i) { return  data[i]; }
    const Element& operator[](const std::size_t i) const { return data[i]; }

    auto begin() {return data.begin();}
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
    auto cbegin() const { return data.cbegin(); }
    auto cend() const { return data.cend(); }

    void push_back(const unsigned index, const float value) { data.emplace_back(Element(index, value)); }
    void push_back(const Element &element) { data.push_back(element); }

    VectorSparse operator*(const float multiplier) const {
        VectorSparse result(*this);
        for (auto& element : result.data) {
            element.value *= multiplier;
        }
        return result;
    }

    VectorSparse operator/(const float divisor) const {
        VectorSparse result(*this);
        for (auto& element : result.data) {
            element.value /= divisor;
        }
        return result;
    }

    VectorSparse operator-() const {
        VectorSparse result(*this);
        for (auto &element : result) {
            element.value = -element.value;
        }
        return result;
    }

    friend VectorSparse operator - (const VectorSparse& a, const VectorSparse& b);
};
