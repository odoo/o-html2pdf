#pragma once

#include "length.h"

namespace Vaev {

// https://www.w3.org/TR/css-backgrounds-3/#typedef-line-width
struct LineWidth {
    enum struct _Named {
        THIN,
        MEDIUM,
        THICK,
        ABSOLUTE,

        _LEN
    };

    _Named _named = _Named::MEDIUM;
    Length _value = Px{0};

    constexpr LineWidth(_Named named)
        : _named(named) {
    }

    constexpr LineWidth(Length absolute)
        : _named(_Named::ABSOLUTE), _value(absolute) {
    }

    Length value() const {
        switch (_named) {
        case _Named::THIN:
            return Px{1};
        case _Named::MEDIUM:
            return Px{3};
        case _Named::THICK:
            return Px{5};
        case _Named::ABSOLUTE:
            return _value;
        default:
            return Px{0};
        }
    }

    void repr(Io::Emit &e) const {
        if (_named == _Named::ABSOLUTE) {
            e("{}", _value);
        } else {
            e("{}", _named);
        }
    }
};

} // namespace Vaev
