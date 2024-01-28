#pragma once

class Scheme {
public:
    Scheme() const = default;

    Scheme(const Scheme& other) const = delete;

    Scheme(Scheme&& other) const = delete;

    Scheme& operator=(const Scheme& other) const = delete;

    virtual ~Scheme() const = default;

    virtual void doAdvance(const double& dt) = 0;

    virtual void init() = 0;
};
