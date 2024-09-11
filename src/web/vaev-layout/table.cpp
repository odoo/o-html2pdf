#include "table.h"

#include "frag.h"

namespace Vaev::Layout {

struct TableFormatingContext {
    struct Col {
        Px width;

        void repr(Io::Emit &e) const {
            e("(column {})", width);
        }
    };

    struct Row {
        Px height;

        void repr(Io::Emit &e) const {
            e("(row {})", height);
        }
    };

    // MARK: Columns -----------------------------------------------------------

    Vec<Col> _cols;
    Px _width;

    void _collectColumnsSizeInner(Tree &t, Frag &f, Input input, usize col) {
        auto display = f.style->display;
        if (display == Display::TABLE_CELL) {
            Input childInput = {
                .commit = Commit::NO,
                .containingBlock = input.containingBlock
            };

            auto childOutput = layout(t, f, childInput);

            if (_cols.len() <= col) {
                _cols.pushBack({
                    .width = childOutput.size.width,
                });
            } else {
                _cols[col].width = max(
                    _cols[col].width,
                    childOutput.size.width
                );
            }
            return;
        }

        for (auto &c : f.children()) {
            Input childInput = {
                .commit = Commit::NO,
                .containingBlock = input.containingBlock,
            };

            if (display.isTableGroup() or display == Display::TABLE) {
                _collectColumnsSizeInner(t, c, childInput, 0);
                continue;
            }

            if (display == Display::TABLE_ROW) {
                _collectColumnsSizeInner(t, c, childInput, col++);
                continue;
            }
        }
    }

    void _collectColumnsSize(Tree &t, Frag &f, Input input) {
        _cols.clear();
        _collectColumnsSizeInner(t, f, input, 0);

        _width = Px{0};
        for (auto &c : _cols)
            _width += c.width;
    }

    // MARK: Rows --------------------------------------------------------------

    Vec<Row> _rows;
    Px _height;

    Px _collectRowSize(Tree &t, Frag &f, Input input) {
        auto display = f.style->display;

        if (display == Display::TABLE_CELL) {
            Input childInput = {
                .commit = Commit::NO,
                .containingBlock = input.containingBlock,
            };
            return layout(t, f, childInput).size.height;
        }

        Px res = Px{0};
        for (auto &c : f.children()) {
            Input childInput = {
                .commit = Commit::NO,
                .availableSpace = Px{0},
                .containingBlock = input.containingBlock,
            };

            res = max(res, _collectRowSize(t, c, childInput));
        }

        return res;
    }

    void _collectRowsSizeInner(Tree &t, Frag &f, Input input) {
        if (not f.style->display.isTableGroup()) {
            _rows.pushBack(Row{_collectRowSize(t, f, input)});
            return;
        }

        for (auto &c : f.children()) {
            Input childInput = {
                .commit = Commit::NO,
                .availableSpace = Px{0},
                .containingBlock = input.containingBlock,
            };

            _collectRowsSizeInner(t, c, childInput);
        }
    }

    void _collectRowsSize(Tree &t, Frag &f, Input input) {
        _rows.clear();
        _collectRowsSizeInner(t, f, input);

        _height = Px{0};
        for (auto &r : _rows)
            _height += r.height;
    }

    // MARK: Layout ------------------------------------------------------------

    void _layoutRow(Tree &t, Frag &f, Input input, usize row) {
        Px inlineSize = {};
        Px rowBlockSize = _rows[row].height;

        for (usize i = 0; i < f.children().len(); i++) {
            auto &c = f.children()[i];

            Input childInput = {
                .commit = input.commit,
                .knownSize = {
                    _cols[i].width,
                    rowBlockSize,
                },
                .availableSpace = {},
                .containingBlock = input.containingBlock,
            };

            auto output = layout(
                t, c,
                childInput
            );

            if (input.commit == Commit::YES) {
                c.layout.position = {
                    inlineSize,
                    Px{0},
                };
            }

            inlineSize += _cols[i].width;
        }
    }

    void _layoutGroup(Tree &t, Frag &f, Input input, usize row) {
        Px blockSize = {};
        Px knownInlineSize = input.knownSize.x.unwrapOr(Px{0});

        for (auto &c : f.children()) {
            auto display = c.style->display;

            Px childBlockSize = {};

            Input childInput = {
                .commit = input.commit,
                .knownSize = {knownInlineSize, NONE},
                .availableSpace = {},
                .containingBlock = {
                    knownInlineSize,
                    Px{0},
                }
            };

            if (display == Display::TABLE_ROW)
                childBlockSize = _collectRowHeight(t, c, childInput);
            else
                childBlockSize = _collectRowsHeight(t, c, childInput);

            childInput.knownSize.height = childBlockSize;

            if (display == Display::TABLE_ROW) {
                _rowLayout(t, c, input, cols);
            } else if (display == Display::TABLE_CAPTION) {
                layout(t, c, childInput);
            } else if (display.isTableGroup()) {
                _placeRows(t, c, input, cols);
            } else {
                layout(t, c, childInput);
            }

            if (input.commit == Commit::YES) {
                c.layout.position = {
                    Px{0},
                    blockSize,
                };
            }

            blockSize += childBlockSize;
        }

        auto inlineSize =
            iter(cols)
                .map([](auto &col) {
                    return col.width;
                })
                .sum();

        return Output::fromSize({
            inlineSize,
            blockSize,
        });
    }
};



Output tableLayout(Tree &t, Frag &f, Input input) {
    TableFormatingContext ctx;
    ctx._collectColumnsSize(t, f, input);
    ctx._collectRowsSize(t, f, input);

    return _placeRows(t, f, input);
}

} // namespace Vaev::Layout
