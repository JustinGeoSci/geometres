// src/counters.rs

//! Unstable diagnostics. Outside semver: may change or disappear in any
//! release.
//!
//! Deterministic operation counters for performance tests. Compiled only
//! under the `unstable-counters` feature. Without it this module does not
//! exist, and every `count!` site in the crate expands to an empty block, so
//! a normal build carries no counter, no thread-local and no branch.
//!
//! # The mechanism
//!
//! One `Cell<u64>` per counter, in a `thread_local!`, incremented by the
//! crate's internal `count!` macro and **read by nothing in the library**: an
//! operation that consulted its own counters would behave differently
//! depending on what ran before it.
//!
//! **A count covers only work done on the calling thread.** [`snapshot`] and
//! [`reset`] read and zero the calling thread's counters and no other's. No
//! operation in this crate spawns a thread, so an operation's whole count
//! lands on the thread that called it, and two tests running in parallel on
//! two threads cannot see each other's counts. A test that reads them resets,
//! runs one operation and snapshots, on one thread; it needs no lock and no
//! `--test-threads=1`.
//!
//! `arithmetica`'s adaptive-ladder counters are a separate mechanism in that
//! crate. A test that reads them depends on `arithmetica` directly, with its
//! `stage-counters` feature; they are process-global, so such a test still
//! runs alone.
//!
//! # The counters
//!
//! | Counter | Counts |
//! |---|---|
//! | `ring_edge_classifications` | ring edges of one areal operand classified against the other in the DE-9IM areal analysis |
//! | `prepared_early_outs` | prepared `relate` calls decided by the envelope early-out |
//! | `prepared_delegations` | prepared `relate` calls that ran the unprepared operation |
//! | `delaunay_orient_calls` | `orient2d` calls made by `ops::delaunay` (walks, ghost tests, collinearity, output certification) |
//! | `delaunay_incircle_calls` | `incircle` calls made by `ops::delaunay` (cavity tests and the canonicalization pass) |
//! | `coverage_simplify_veto_checks` | removals `ops::coverage_simplify` put to its exact topology checks |
//! | `coverage_simplify_veto_rejections` | of those, the removals the checks rejected |
//! | `coverage_simplify_index_updates` | segment insertions and removals in `ops::coverage_simplify`'s dynamic index |
//! | `coverage_simplify_index_queries` | candidate queries `ops::coverage_simplify` made of that index |

use core::cell::Cell;

/// Defines the counter set once: the thread-local block, its snapshot type,
/// and the read and reset functions, so a counter cannot be added to one and
/// forgotten in another.
macro_rules! define_counters {
    ($($(#[doc = $doc:literal])* $name:ident),* $(,)?) => {
        /// One thread's counters. One `Cell<u64>` per counter.
        pub(crate) struct Counters {
            $(pub(crate) $name: Cell<u64>,)*
        }

        thread_local! {
            /// The calling thread's instance. `const`-initialized and without
            /// a destructor, so it is never in a destroyed state while its
            /// thread runs.
            pub(crate) static COUNTERS: Counters = const {
                Counters { $($name: Cell::new(0),)* }
            };
        }

        /// A copy of every counter of the calling thread at one moment.
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub struct CounterSnapshot {
            $($(#[doc = $doc])* pub $name: u64,)*
        }

        /// Reads every counter of the calling thread.
        #[must_use]
        pub fn snapshot() -> CounterSnapshot {
            COUNTERS
                .try_with(|c| CounterSnapshot {
                    $($name: c.$name.get(),)*
                })
                .unwrap_or_default()
        }

        /// Zeroes every counter of the calling thread.
        pub fn reset() {
            let _ = COUNTERS.try_with(|c| {
                $(c.$name.set(0);)*
            });
        }
    };
}

define_counters! {
    /// Ring edges classified by `summarize_segment` in the areal DE-9IM
    /// analysis (`polygon_analysis::ring_summary`): every ring edge of A
    /// against B, then every ring edge of B against A, when the operands'
    /// boxes meet. `O(n·m)` work in total; prepared `relate` lowers it.
    ring_edge_classifications,
    /// Prepared `relate` calls answered by the envelope early-out alone.
    prepared_early_outs,
    /// Prepared `relate` calls that ran the unprepared operation: a
    /// delegating (non-simple) mask, or an empty operand.
    prepared_delegations,
    /// `orient2d` calls made by `ops::delaunay`: point location walks,
    /// ghost-triangle tests, collinearity of the starting chain, and the
    /// output certification.
    delaunay_orient_calls,
    /// `incircle` calls made by `ops::delaunay`: cavity membership during
    /// construction, and the canonicalization pass over internal edges.
    delaunay_incircle_calls,
    /// Removals `ops::coverage_simplify` put to its exact topology checks
    /// (the veto): every candidate vertex that is not a node.
    coverage_simplify_veto_checks,
    /// Removals the topology checks rejected; the vertex is then kept.
    coverage_simplify_veto_rejections,
    /// Segment insertions and removals in `ops::coverage_simplify`'s dynamic
    /// segment index: the initial load, and three per accepted removal.
    coverage_simplify_index_updates,
    /// Candidate queries `ops::coverage_simplify` made of its dynamic
    /// segment index: one per topology check that reaches it.
    coverage_simplify_index_queries,
}

/// Adds `n` to one counter of the calling thread. Crate-internal; used only
/// by the `count!` macro. Wrapping, as the atomic counters it replaces were.
#[inline]
pub(crate) fn add(field: fn(&Counters) -> &Cell<u64>, n: u64) {
    let _ = COUNTERS.try_with(|c| {
        let cell = field(c);
        cell.set(cell.get().wrapping_add(n));
    });
}
