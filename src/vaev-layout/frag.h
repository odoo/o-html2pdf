#pragma once

#include <karm-image/picture.h>
#include <karm-text/run.h>
#include <vaev-markup/dom.h>
#include <vaev-style/computer.h>

#include "base.h"

namespace Vaev::Layout {

// MARK: Frag ------------------------------------------------------------------

using Content = Union<
    None,
    Vec<Frag>,
    Strong<Text::Run>,
    Image::Picture>;

struct Frag : public Meta::NoCopy {
    Strong<Style::Computed> style;
    Text::Font font;
    Content content = NONE;
    Layout layout;
    Cow<TableSpan> tableSpan;

    Frag(Strong<Style::Computed> style, Text::Font font);

    Frag(Strong<Style::Computed> style, Text::Font font, Content content);

    Slice<Frag> children() const;

    MutSlice<Frag> children();

    void add(Frag &&frag);

    void repr(Io::Emit &e) const;
};

struct Tree {
    Frag root;
    Viewport viewport;
};

// MARK: Build -----------------------------------------------------------------

void _buildNode(Style::Computer &c, Markup::Node const &n, Frag &parent);

Frag build(Style::Computer &c, Markup::Document const &doc);

// MARK: Layout ----------------------------------------------------------------

InsetsPx computeMargins(Tree &t, Frag &f, Input input);

Output layout(Tree &t, Frag &f, Input input);

Px measure(Tree &t, Frag &f, Axis axis, IntrinsicSize intrinsic, Px availableSpace);

void wireframe(Frag &frag, Gfx::Canvas &g);

} // namespace Vaev::Layout
