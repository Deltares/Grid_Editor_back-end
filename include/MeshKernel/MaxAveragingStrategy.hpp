﻿#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief MaxAveragingStrategy implements the AveragingStrategy which takes the maximum value of the added points.
    class MaxAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new MaxAveragingStrategy.
        /// @param[in] missingValue The value used to indicate a missing value.
        explicit MaxAveragingStrategy(double missingValue);

        void Add(Point const& samplePoint, double sampleValue) override;

        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The current result returned in Calculate.
        double m_result = std::numeric_limits<double>::lowest();

        /// @brief The value returned when no valid value can be returned.
        double const m_missingValue;
    };
} // namespace meshkernel::averaging
