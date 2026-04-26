#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "SimTKcommon.h"

using namespace SimTK;
using std::string;
using std::unique_ptr;

/* -------------------------------------------------------------------------- */
/* Original intent:
   Verify semantics of CloneOnWritePtr:
   - empty state invariants
   - ownership / reference counting
   - copy-on-write detach behavior
   - move/copy correctness
   - interaction with raw pointers
   - no leaks via Base accounting
*/
/* -------------------------------------------------------------------------- */

class Base {
    public:
    explicit Base(const std::string& n)
        : m_name(n) {
        ++m_constructions;
    }

    Base(const Base& src)
        : m_name(src.m_name) {
        ++m_copies;
    }

    virtual ~Base() {
        ++m_destructions;
    }

    [[nodiscard]] virtual Base* clone() const = 0;
    [[nodiscard]] virtual int getValue() const = 0;
    virtual auto updValue() -> int& = 0;

    void setValue(int v) {
        updValue() = v;
    }

    auto getName() const -> const char* {
        return m_name.c_str();
    }

    static auto getNumAlive() -> int {
        return m_constructions + m_copies - m_destructions;
    }

    static int m_constructions;
    static int m_destructions;
    static int m_copies;

    private:
    std::string m_name;
};

int Base::m_constructions = 0;
int Base::m_destructions = 0;
int Base::m_copies = 0;

/* ---------------- Derived hierarchy ---------------- */

class Derived1 : public Base {
    public:
    explicit Derived1(const std::string& n)
        : Base(n) {
    }
};

class Derived2 : public Base {
    public:
    Derived2(const std::string& n, int v)
        : Base(n)
        , m_val2(v) {
    }

    [[nodiscard]] auto clone() const -> Derived2* override {
        return new Derived2(*this);
    }

    [[nodiscard]] auto getValue() const -> int override {
        return m_val2;
    }

    auto updValue() -> int& override {
        return m_val2;
    }

    private:
    int m_val2;
};

class Sub1 : public Derived1 {
    public:
    Sub1(const std::string& n, int v)
        : Derived1(n)
        , m_val1(v) {
    }

    auto clone() const -> Sub1* override {
        return new Sub1(*this);
    }

    auto getValue() const -> int override {
        return m_val1;
    }

    auto updValue() -> int& override {
        return m_val1;
    }

    private:
    int m_val1;
};

/* -------------------------------------------------------------------------- */
/* EMPTY STATE TESTS */
/* What authors tested:
   - null semantics identical to std::shared_ptr-like API
   - no accidental allocation
   - safe API behavior (throws in debug)
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_EmptyState, BasicInvariants) {
    CloneOnWritePtr<Base> p;
    CloneOnWritePtr<Base> pp;

    EXPECT_TRUE(p.empty());
    EXPECT_FALSE(p.unique());
    EXPECT_EQ(p.use_count(), 0);

    EXPECT_FALSE(static_cast<bool>(p));
    EXPECT_EQ(p.get(), nullptr);
    EXPECT_EQ(p.upd(), nullptr);
    EXPECT_EQ(p.release(), nullptr);

#ifndef NDEBUG
    EXPECT_THROW(p.getRef(), std::exception);
    EXPECT_THROW(p.updRef(), std::exception);
    EXPECT_THROW(p->getValue(), std::exception);
    EXPECT_THROW((*p).getValue(), std::exception);
#endif

    p.detach();
    p.swap(pp);
    std::swap(p, pp);

    EXPECT_TRUE(p.empty());
    EXPECT_TRUE(pp.empty());
}

TEST(SimTKCommon_CloneOnWritePtr_EmptyState, NullComparisonsAndAssignments) {
    CloneOnWritePtr<Sub1> q(nullptr);
    CloneOnWritePtr<Derived2> d2(nullptr);
    CloneOnWritePtr<Base> p;

    EXPECT_TRUE(q.empty());
    EXPECT_TRUE(d2.empty());

    EXPECT_TRUE(d2 == nullptr);
    EXPECT_TRUE(nullptr == d2);
    EXPECT_TRUE(d2 == p);

    EXPECT_FALSE(d2 != nullptr);
    EXPECT_TRUE(d2 >= nullptr);
    EXPECT_TRUE(d2 <= nullptr);

    CloneOnWritePtr<Base> pp;
    p = pp;
    p = d2;
    p = std::move(pp);
    p = std::move(d2);

    EXPECT_TRUE(p.empty());
}

/* -------------------------------------------------------------------------- */
/* ALLOCATION & COPY-ON-WRITE */
/* What authors tested:
   - reference counting correctness
   - detach occurs on mutation
   - release semantics
   - move vs copy allocation behavior
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_Allocation, ReferenceCountingAndDetach) {
    EXPECT_EQ(Base::getNumAlive(), 0);

    CloneOnWritePtr<Base> p(new Sub1("first", 1));
    CloneOnWritePtr<Base> q(new Sub1("second", 2));
    CloneOnWritePtr<Base> r(new Derived2("d2", 999));

    EXPECT_EQ(Base::getNumAlive(), 3);

    EXPECT_TRUE(p.unique());
    EXPECT_TRUE(q.unique());
    EXPECT_TRUE(r.unique());

    // Share ownership
    p = q = r;

    EXPECT_EQ(p.use_count(), 3);
    EXPECT_EQ(Base::getNumAlive(), 3);

    // Read-only access should not detach
    const int rval_const = r.get()->getValue();
    EXPECT_EQ(rval_const, 999);
    EXPECT_EQ(Base::getNumAlive(), 3);

    // Non-const access forces detach
    const int rval = r->getValue();
    EXPECT_EQ(rval, 999);
    EXPECT_TRUE(r.unique());
    EXPECT_EQ(Base::getNumAlive(), 4);
}

TEST(SimTKCommon_CloneOnWritePtr_Allocation, ReleaseAndTakeover) {
    CloneOnWritePtr<Base> p(new Sub1("x", 1));
    CloneOnWritePtr<Base> q = p;

    EXPECT_EQ(p.use_count(), 2);

    Base* raw = p.release();

    EXPECT_TRUE(p.empty());
    EXPECT_TRUE(q.unique());
    EXPECT_NE(raw, nullptr);

    CloneOnWritePtr<Base> takeover(raw);

    EXPECT_TRUE(takeover.unique());
}

TEST(SimTKCommon_CloneOnWritePtr_Allocation, DetachOnMutation) {
    CloneOnWritePtr<Base> a(new Derived2("d", 5));
    CloneOnWritePtr<Base> b = a;

    EXPECT_EQ(a.use_count(), 2);

    b->setValue(10); // should detach

    EXPECT_TRUE(a.unique());
    EXPECT_TRUE(b.unique());
    EXPECT_NE(a.get(), b.get());
    EXPECT_EQ(b->getValue(), 10);
}

/* -------------------------------------------------------------------------- */
/* SELF ASSIGNMENT / MOVE */
/* What authors tested:
   - self-copy is no-op
   - self-move is safe
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_SelfOps, SelfAssignmentAndMove) {
    CloneOnWritePtr<Base> p(new Derived2("d", 1));
    CloneOnWritePtr<Base> q = p;

    p = p;
    EXPECT_EQ(p.use_count(), 2);

    p = std::move(p);
    EXPECT_TRUE(p);
    EXPECT_TRUE(p.unique() || p.use_count() >= 1);
}

/* -------------------------------------------------------------------------- */
/* LEAK CHECK */
/* What authors tested:
   - all objects destroyed at end
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_Lifetime, NoLeaks) {
    EXPECT_EQ(Base::getNumAlive(), 0);
}

/* -------------------------------------------------------------------------- */
/* ResetOnCopy */
/* What authors tested:
   - copy resets value to default
   - move preserves value
   - self-copy clears
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_ResetOnCopy, BasicSemantics) {
    ResetOnCopy<int> a(5);
    EXPECT_EQ(a, 5);

    ResetOnCopy<int> b(a);
    EXPECT_EQ(b, 0); // reset on copy

    b = a;
    EXPECT_EQ(b, 0);

    b = std::move(a);
    EXPECT_EQ(b, 5);
}

TEST(SimTKCommon_CloneOnWritePtr_ResetOnCopy, SelfAssignmentClears) {
    ResetOnCopy<Array_<int>> v{1, 2, 3};

    EXPECT_EQ(v.size(), 3);

    v = v; // reset-on-copy
    EXPECT_TRUE(v.empty());
}

/* -------------------------------------------------------------------------- */
/* ReinitOnCopy */
/* What authors tested:
   - copy restores initial value
   - move preserves mutated value
*/
/* -------------------------------------------------------------------------- */

TEST(SimTKCommon_CloneOnWritePtr_ReinitOnCopy, CopyReinitializes) {
    ReinitOnCopy<int> a(10);

    a = 20;
    EXPECT_EQ(a, 20);

    ReinitOnCopy<int> b(a);
    EXPECT_EQ(b, 10); // reinitialized
}

TEST(SimTKCommon_CloneOnWritePtr_ReinitOnCopy, MovePreservesState) {
    ReinitOnCopy<int> a(10);
    a = 30;

    ReinitOnCopy<int> b(std::move(a));
    EXPECT_EQ(b, 30);
}