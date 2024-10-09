#include "table.h"

#include "frag.h"
#include "values.h"

namespace Vaev::Layout {

bool isHeadBodyFootOrRow(Display display) {
    return (
        display == Vaev::Display::Internal::TABLE_HEADER_GROUP or
        display == Vaev::Display::Internal::TABLE_ROW_GROUP or
        display == Vaev::Display::Internal::TABLE_FOOTER_GROUP or
        display == Vaev::Display::Internal::TABLE_ROW
    );
}

bool isCaseOfStep8(Display display) {
    return (
        isHeadBodyFootOrRow(display) or
        display == Vaev::Display::Internal::TABLE_COLUMN_GROUP
    );
}

Frag *nextInRange(Frag *curr, Frag *end, Func<bool(Display)> pred) {
    while (curr != end) {
        if (pred(curr->style->display))
            return curr;
        curr++;
    }
    return curr;
}

struct Slots {
    struct Cell {
        Math::Vec2u anchor;
        Frag *el;
    };

    Vec<Vec<Cell *>> rows = {{nullptr}};
    Math::Vec2u size = {1, 1};

    void increaseWidth(usize span = 1) {

        for (auto &row : rows) {
            for (usize i = 0; i < span; ++i)
                row.pushBack(nullptr);
        }
        size.x += span;
    }

    void increaseHeight(usize span = 1) {

        Vec<Cell *> newRow{Buf<Cell *>::init(size.x, nullptr)};
        for (usize i = 0; i < span; ++i) {
            rows.pushBack(newRow);
        }
        size.y += span;
    }

    Cell *get(usize x, usize y) {
        // FIXME: maybe remove in prod? guaranteed that algo is correct
        if (x >= size.x or y >= size.y)
            panic("bad coordinates for table slot");

        return rows[y][x];
    }

    void set(usize x, usize y, Cell *cell) {

        // FIXME: maybe remove in prod? guaranteed that algo is correct
        if (x >= size.x or y >= size.y)
            panic("bad coordinates for table slot");

        rows[y][x] = cell;
    }
};

struct TableForming {
    Math::Vec2u current;

    struct ColCorresp {
        usize start, end;
        Frag const *el;
    };

    struct Group {
        enum {
            ROW,
            COLUMN
        } type;

        usize start, end;
        Frag const *el;
    };

    Vec<ColCorresp> colCorresps;
    Vec<Group> groups;

    Slots &slots;

    TableForming(Slots &s) : slots(s) {}

    struct DownwardsGrowingCell {
        Slots::Cell *cell;
        usize xpos, width;
    };

    Vec<DownwardsGrowingCell> downwardsGrowingCells;

    void processRows(Frag *tableRowElement) {
        auto tableRowChildren = tableRowElement->children();

        if (slots.size.y == current.y)
            slots.increaseHeight();

        current.x = 0;

        // Run the algorithm for growing downward-growing cells.

        auto currTdOrThChild = nextInRange(tableRowChildren.begin(), tableRowChildren.end(), [](Display d) {
            return d == Display::TABLE_CELL;
        });

        while (currTdOrThChild != tableRowChildren.end()) {

            while (current.x < slots.size.x and slots.get(current.x, current.y))
                current.x++;

            if (current.x == slots.size.x)
                slots.increaseWidth();

            usize rowSpan = currTdOrThChild->tableSpan->row, colSpan = currTdOrThChild->tableSpan->col;
            bool cellGrowsDownward;
            if (rowSpan == 0 /* and the table element's node document is not set to quirks mode, */) {
                cellGrowsDownward = true;
                rowSpan = 1;
            } else
                cellGrowsDownward = false;

            if (slots.size.x < current.x + colSpan) {
                slots.increaseWidth(current.x + colSpan - slots.size.x);
            }

            if (slots.size.y < current.y + rowSpan) {
                slots.increaseHeight(current.y + rowSpan - slots.size.y);
            }

            // 13
            {
                Slots::Cell *cell = new Slots::Cell{
                    .anchor = current,
                    .el = currTdOrThChild
                };

                for (usize x = current.x; x < current.x + colSpan; ++x) {
                    for (usize y = current.y; y < current.y + rowSpan; ++y) {
                        slots.set(x, y, cell);
                    }
                }

                if (cellGrowsDownward) {
                    downwardsGrowingCells.pushBack({.cell = cell, .xpos = current.x, .width = colSpan});
                }
            }

            current.x += colSpan;

            currTdOrThChild++;
            currTdOrThChild = nextInRange(currTdOrThChild, tableRowChildren.end(), [](Display d) {
                return d == Display::TABLE_ROW or d == Display::TABLE_CELL;
            });
        }
        current.y++;
    }

    void growDownwardGrowingCells() {
        for (auto &[cell, cellx, width] : downwardsGrowingCells) {
            for (usize x = cellx; x < cellx + width; x++)
                slots.set(x, current.y, cell);
        }
    }

    void endRowGroup() {
        while (current.y < slots.size.y) {
            growDownwardGrowingCells();
            current.y++;
            downwardsGrowingCells.clear();
        }
    }

    void processRowGroup(Frag *rowGroupElement) {
        usize ystart = slots.size.y;

        auto nextTrElement = nextInRange(rowGroupElement->children().begin(), rowGroupElement->children().end(), [](Display d) {
            return d == Display::TABLE_ROW;
        });
        while (nextTrElement != rowGroupElement->children().end()) {

            processRows(nextTrElement);

            nextTrElement++;
            nextTrElement = nextInRange(nextTrElement, rowGroupElement->children().end(), [](Display d) {
                return d == Display::TABLE_ROW;
            });
        }

        if (slots.size.y > ystart) {
            groups.pushBack({.type = Group::ROW, .start = ystart, .end = slots.size.y - 1, .el = rowGroupElement});
        }

        endRowGroup();
    }

    // https://html.spec.whatwg.org/multipage/tables.html#forming-a-table
    void run(Frag *tableElement) {
        auto currTableChild = tableElement->children().begin();
        auto endTableChild = tableElement->children().end();

        if (currTableChild == endTableChild)
            return;

        // TODO: Associate the first caption element child of the table element with the table. If there are no such children, then it has no associated caption element.

        currTableChild = nextInRange(currTableChild, endTableChild, [](Display d) {
            return isCaseOfStep8(d);
        });

        if (currTableChild == endTableChild)
            return;

        // MARK: Columns groups
        while (currTableChild->style->display == Display::TABLE_COLUMN_GROUP) {

            auto colEl = nextInRange(currTableChild->children().begin(), currTableChild->children().end(), [](Display d) {
                return d == Vaev::Display::TABLE_COLUMN;
            });
            if (colEl) {
                usize startColRange = slots.size.x;

                // MARK: Columns
                while (colEl) {
                    auto span = colEl->tableSpan->col;
                    slots.increaseWidth(span);

                    colCorresps.pushBack({.start = slots.size.x - span, .end = slots.size.x - 1, .el = colEl});

                    colEl = nextInRange(colEl + 1, currTableChild->children().end(), [](Display d) {
                        return isCaseOfStep8(d);
                    });
                }

                groups.pushBack({.type = Group::COLUMN, .start = startColRange, .end = slots.size.x - 1, .el = currTableChild});
            } else {
                auto span = currTableChild->tableSpan->col;
                slots.increaseWidth(span);

                groups.pushBack({.type = Group::COLUMN, .start = slots.size.x - span + 1, .end = slots.size.x - 1, .el = currTableChild});
            }

            currTableChild++;
            currTableChild = nextInRange(currTableChild, endTableChild, [](Display d) {
                return isCaseOfStep8(d);
            });
        }

        current.y = 0;

        // MARK: rows
        while (true) {

            currTableChild = nextInRange(currTableChild, endTableChild, [](Display d) {
                return isHeadBodyFootOrRow(d);
            });

            if (currTableChild == endTableChild)
                break;

            if (currTableChild->style->display == Display::TABLE_ROW) {
                processRows(currTableChild);
                currTableChild++;
                continue;
            }

            endRowGroup();

            if (currTableChild->style->display == Display::TABLE_FOOTER_GROUP) {
                // If the current element is a tfoot, then add that element to the list of pending tfoot elements,
                // advance the current element to the next child of the table, and return to the step labeled rows.
                currTableChild++;
                continue;
            }

            // The current element is either a thead or a tbody.
            if (currTableChild->style->display != Display::TABLE_HEADER_GROUP and currTableChild->style->display != Display::TABLE_ROW_GROUP) {
                // FIXME: prod code should not fail, but ok for current dev scenario
                panic("current element should be thead or tbody");
            }

            processRowGroup(currTableChild);

            currTableChild++;
        }
    }
};

Tuple<Vec<Px>, Px> getFixedColWidths(Tree &t, Frag &f, Input &input, Slots &s, TableForming &table) {
    // FIXME: should be considering COL elements and not COL GROUPS?
    Vec<Px> colWidth{Buf<Px>::init(s.size.x, Px{0})};
    for (auto &group : table.groups) {
        if (group.type != TableForming::Group::COLUMN)
            continue;
        auto width = group.el->style->sizing->width;
        if (width == Size::Type::AUTO)
            continue;

        for (usize x = group.start; x <= group.end; ++x) {
            colWidth[x] = resolve(t, f, width.value, input.availableSpace.x);
        }
    }

    // Using first row cells to define columns widths

    usize x = 0;
    while (x < s.size.x) {
        if (colWidth[x] == Px{0}) {
            auto *cell = s.get(x, 1);

            if (cell->el->style->sizing->width == Size::Type::AUTO) {
                x++;
                continue;
            }

            auto cellWidth = resolve(t, f, cell->el->style->sizing->width.value, input.availableSpace.x);
            auto colSpan = cell->el->tableSpan->col;

            for (usize j = 0; j < colSpan; ++j, x++) {
                // FIXME: not overriding values already computed, but should we subtract the already computed from
                // cellWidth before division?
                if (colWidth[x] == Px{0})
                    colWidth[x] = cellWidth / Px{colSpan};
            }
        }
        x++;
    }

    // FIXME:  (plus cell spacing or borders)
    Px sumColsWidths{0};
    usize emptyCols{0};
    for (auto w : colWidth) {
        if (w == Px{0})
            emptyCols++;
        else
            sumColsWidths += w;
    }

    if (sumColsWidths < input.availableSpace.x) {
        Px toDistribute = sumColsWidths / Px{emptyCols};
        for (auto &w : colWidth)
            if (w == Px{0})
                w = toDistribute;
        sumColsWidths = input.availableSpace.x;
    }

    // This method only runs in the FIXED table layout algorithm, where we should have a definite table width
    if (resolve(t, f, f.style->sizing->width.value, Px{0})) {
        Px toDistribute = input.knownSize.x.unwrap() / Px{s.size.x};
        for (auto &w : colWidth)
            w += toDistribute;
        sumColsWidths = input.knownSize.x.unwrap();
    }
    // MARK: Fixed Table Layout

    return {colWidth, sumColsWidths};
}

// MARK: Auto Table Layout
Tuple<Vec<Px>, Px> getAutoColWidths(Tree &t, Frag &f, Input input, Slots &s, TableForming &table) {
    // https://www.w3.org/TR/CSS21/tables.html#auto-table-layout

    Vec<Px> minColWidth{Buf<Px>::init(s.size.x, Px{0})};

    Vec<Px> maxColWidth{Buf<Px>::init(s.size.x, Px{0})};

    for (usize i = 0; i < s.size.y; ++i) {
        for (usize j = 0; j < s.size.x; ++j) {
            auto cell = s.get(j, i);

            // FIXME
            if (cell == nullptr)
                continue;

            if (not(cell->anchor == Math::Vec2u{j, i}))
                continue;

            auto colSpan = cell->el->tableSpan->col;
            if (colSpan > 1)
                continue;

            // FIXME: what should be the parameter for intrinsic in the vertical axis?
            auto cellMinOutput = layout(
                t,
                *cell->el,
                Input{
                    .commit = Commit::NO,
                    .intrinsic = {IntrinsicSize::MIN_CONTENT, IntrinsicSize::AUTO},
                    .knownSize = {NONE, NONE}
                }
            );

            auto cellMinWidth = cellMinOutput.size.x;
            if (cell->el->style->sizing->width != Size::Type::AUTO) {
                cellMinWidth = max(
                    cellMinWidth,
                    resolve(t, f, cell->el->style->sizing->width.value, Px{0})
                );
            }

            minColWidth[j] = max(minColWidth[j], cellMinWidth);

            // FIXME: what should be the parameter for intrinsic in the vertical axis?
            auto cellMaxOutput = layout(
                t,
                *cell->el,
                Input{
                    .commit = Commit::NO,
                    .intrinsic = {IntrinsicSize::MAX_CONTENT, IntrinsicSize::AUTO},
                    .knownSize = {NONE, NONE}
                }
            );

            maxColWidth[j] = max(maxColWidth[j], cellMaxOutput.size.x);
        }
    }

    for (auto &[start, end, el] : table.colCorresps) {
        auto width = el->style->sizing->width;
        if (width == Size::Type::AUTO)
            continue;
        auto widthValue = resolve(t, f, width.value, input.availableSpace.x);

        for (usize x = start; x <= end; ++x) {
            minColWidth[x] = max(minColWidth[x], widthValue);
            maxColWidth[x] = max(maxColWidth[x], widthValue);
        }
    }

    for (usize i = 0; i < s.size.y; ++i) {
        for (usize j = 0; j < s.size.x; ++j) {
            auto cell = s.get(j, i);

            // FIXME
            if (cell == nullptr)
                continue;

            if (not(cell->anchor == Math::Vec2u{j, i}))
                continue;

            auto colSpan = cell->el->tableSpan->col;
            if (colSpan <= 1)
                continue;

            {
                // FIXME: what should be the parameter for intrinsic in the vertical axis?
                auto cellMinOutput = layout(
                    t,
                    *cell->el,
                    Input{
                        .commit = Commit::NO,
                        .intrinsic = {IntrinsicSize::MIN_CONTENT, IntrinsicSize::AUTO},
                        .knownSize = {NONE, NONE}
                    }
                );

                auto cellMinWidthContribution = cellMinOutput.size.x / Px{colSpan};
                // FIXME
                if (cellMinWidthContribution > Px{0}) {
                    for (usize k = 0; k < colSpan; ++k) {

                        if (minColWidth[j + k] == Limits<Px>::MAX)
                            minColWidth[j + k] = cellMinWidthContribution;
                        else
                            minColWidth[j + k] += cellMinWidthContribution;
                    }
                }
            }
            {

                // FIXME: what should be the parameter for intrinsic in the vertical axis?
                auto cellMaxOutput = layout(
                    t,
                    *cell->el,
                    Input{
                        .commit = Commit::NO,
                        .intrinsic = {IntrinsicSize::MAX_CONTENT, IntrinsicSize::AUTO},
                        .knownSize = {NONE, NONE}
                    }
                );
                auto cellMaxWidthContribution = cellMaxOutput.size.x / Px{colSpan};
                for (usize k = 0; k < colSpan; ++k) {
                    if (maxColWidth[j + k] == Limits<Px>::MAX)
                        maxColWidth[j + k] = cellMaxWidthContribution;
                    else
                        maxColWidth[j + k] += cellMaxWidthContribution;
                }
            }
        }
    }

    for (auto &group : table.groups) {
        if (group.type != TableForming::Group::COLUMN)
            continue;

        auto columnGroupWidth = group.el->style->sizing->width;
        if (columnGroupWidth == Size::Type::AUTO)
            continue;

        Px currSumOfGroupWidth{0};
        for (usize x = group.start; x <= group.end; ++x) {
            currSumOfGroupWidth += minColWidth[x];
        }

        auto columnGroupWidthValue = resolve(t, f, columnGroupWidth.value, input.availableSpace.x);
        if (currSumOfGroupWidth >= columnGroupWidthValue)
            continue;

        Px toDistribute = (columnGroupWidthValue - currSumOfGroupWidth) / Px{group.end - group.start + 1};
        for (usize x = group.start; x <= group.end; ++x) {
            minColWidth[x] += toDistribute;
        }
    }

    // FIXME
    for (auto &x : minColWidth)
        x = max(x, Px{15});
    for (auto &x : maxColWidth)
        x = max(x, Px{15});

    Px sumMinColWidths{0}, sumMaxColWidths{0};
    for (auto x : minColWidth)
        sumMinColWidths += x;
    for (auto x : maxColWidth)
        sumMaxColWidths += x;

    // FIXME: implement CAPMIN, no need to refactor now since will change soon
    if (f.style->sizing->width == Size::Type::AUTO) {
        if (sumMaxColWidths < input.containingBlock.x) {
            return {maxColWidth, sumMaxColWidths};
        } else {
            return {minColWidth, sumMinColWidths};
        }
    } else {
        auto tableComputedWidth = resolve(t, f, f.style->sizing->width.value, Px{0});
        auto tableUsedWidth = max(tableComputedWidth, sumMinColWidths);

        if (sumMinColWidths < tableUsedWidth) {
            auto toDistribute = (tableUsedWidth - sumMinColWidths) / Px{s.size.x};
            for (auto &w : minColWidth)
                w += toDistribute;
        }
        return {minColWidth, tableUsedWidth};
    }
}

Output tableLayout(Tree &t, Frag &f, Input input) {

    Slots s;

    TableForming table(s);
    table.run(&f);

    // TODO: its also possible to run the fixed version even if the width is auto, to be discussed
    bool shouldRunAutoAlgorithm = *f.style->tableLayout == TableLayout::AUTO or f.style->sizing->width == Size::AUTO;

    auto [colWidth, sumColsWidths] = shouldRunAutoAlgorithm
                                         ? getAutoColWidths(t, f, input, s, table)
                                         : getFixedColWidths(t, f, input, s, table);

    Vec<Px> rowHeight{Buf<Px>::init(s.size.y, Px{0})};
    { // Using table-column and table-row elements dimensions to compute the rows and cols dimensions
        for (auto &group : table.groups) {
            if (group.type != TableForming::Group::ROW)
                continue;

            auto height = group.el->style->sizing->height;
            if (height == Size::Type::AUTO)
                continue;

            for (usize y = group.start; y <= group.end; ++y) {
                rowHeight[y] = resolve(t, f, height.value, Px{0});
            }
        }
    }

    for (usize i = 1; i < s.size.y; ++i) {
        for (usize j = 0; j < s.size.x; ++j) {
            auto cell = s.get(j, i);

            // FIXME
            if (cell == nullptr)
                continue;

            if (not(cell->anchor == Math::Vec2u{j, i}))
                continue;

            auto rowSpan = cell->el->tableSpan->row;

            if (cell->el->style->sizing->height != Size::Type::AUTO) {
                auto computedHeight = resolve(t, f, cell->el->style->sizing->height.value, Px{0});

                // the computed height of each cell spanning the current row exclusively (if definite, percentages being treated as 0px), and
                for (usize k = 0; k < rowSpan; k++) {
                    rowHeight[i + k] = max(rowHeight[i + k], Px{computedHeight / Px{rowSpan}});
                }
            }

            auto cellOutput = layout(
                t,
                *cell->el,
                Input{
                    .commit = Commit::NO,
                    .intrinsic = {IntrinsicSize::AUTO, IntrinsicSize::MIN_CONTENT},
                    .knownSize = {colWidth[j], NONE}
                }
            );
            // the computed height of each cell spanning the current row exclusively (if definite, percentages being treated as 0px), and
            for (usize k = 0; k < rowSpan; k++) {
                rowHeight[i + k] = max(rowHeight[i + k], Px{cellOutput.size.y / Px{rowSpan}});
            }
        }
    }

    if (input.commit == Commit::YES) {
        auto buildPref = [](Vec<Px> const &v) {
            Vec<Px> pref(v);
            for (usize i = 1; i < v.len(); ++i) {
                pref[i] = pref[i - 1] + pref[i];
            }
            return pref;
        };
        auto queryPref = [](Vec<Px> const &pref, usize l, usize r) {
            if (l == 0)
                return pref[r];
            return pref[r] - pref[l - 1];
        };

        auto colWidthPref = buildPref(colWidth), rowHeightPref = buildPref(rowHeight);

        Px currPositionY{0};
        for (usize i = 1; i < s.size.y; currPositionY += rowHeight[i], i++) {
            Px currPositionX = Px{0};
            for (usize j = 0; j < s.size.x; currPositionX += colWidth[j], j++) {
                auto cell = s.get(j, i);
                // FIXME
                if (cell == nullptr)
                    continue;

                if (not(cell->anchor == Math::Vec2u{j, i}))
                    continue;

                auto colSpan = cell->el->tableSpan->col;
                auto rowSpan = cell->el->tableSpan->row;

                auto cellOutput = layout(
                    t,
                    *cell->el,
                    Input{
                        .commit = Commit::YES,
                        .knownSize = {
                            queryPref(colWidthPref, j, j + colSpan - 1),
                            queryPref(rowHeightPref, i, i + rowSpan - 1)
                        },
                        .position{currPositionX, currPositionY}
                    }
                );
            };
        }
    }

    Px sumRowsHeights{0};
    for (auto h : rowHeight)
        sumRowsHeights += h;

    return Output::fromSize({sumColsWidths, sumRowsHeights});
}

} // namespace Vaev::Layout
