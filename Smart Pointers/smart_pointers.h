#include <iostream>
#include <memory>
#include <optional>


struct Counter {
    size_t shared_counter;
    size_t weak_counter;

    Counter(size_t shared_counter, size_t weak_counter):
            shared_counter(shared_counter), weak_counter(weak_counter) {}
};

struct ControlBlockBase {
    virtual Counter& getCounter() = 0;

    virtual void destroy_data(void*) = 0;
    virtual void destroy_partial() = 0;
    virtual void destroy_all() = 0;
};

template<typename T, typename Deleter = std::default_delete<T>, typename Allocator = std::allocator<T>>
class ControlBlock: public ControlBlockBase {
private:
    std::optional<T> data;
    Counter counter;
    Allocator object_allocator;
    Deleter object_deleter;

    template <typename U>
    friend class SharedPtr;
public:
    bool need_destuction = false;

    using ControlBlockType = ControlBlock<T, Deleter, Allocator>;
    using ControlBlockAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<ControlBlockType>;
    using AllocatorTraits = std::allocator_traits<ControlBlockAllocator>;

    void destroy_data(void* data) override {
        object_deleter(static_cast<T*>(data));
    }

    void destroy_partial() override {
        data.reset();
    }

    ControlBlock(): counter(1, 0) {}

    explicit ControlBlock(const Deleter& destroy_data):
            counter(1, 0),
            object_deleter(destroy_data) {}

    ControlBlock(const Deleter& destroy_data, const Allocator& allocator):
            counter(1, 0),
            object_allocator(allocator),
            object_deleter(destroy_data) {}

    template <typename... Args>
    explicit ControlBlock(const Allocator& allocator, Args&&... args):
            data(std::in_place_t(), std::forward<Args>(args)...),
            counter(1, 0),
            object_allocator(allocator) {}

    Counter& getCounter() override {
        return counter;
    }

    void destroy_all() override {
        ControlBlockAllocator save(object_allocator);
        if (need_destuction)
            AllocatorTraits::destroy(save, this);
        AllocatorTraits::deallocate(save, this, 1);
    }
};

template <typename T>
class WeakPtr;

template <typename T>
class SharedPtr {
private:
    T* data = nullptr;
    ControlBlockBase* cptr = nullptr;

public:
    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

    template <typename U, typename Allocator, typename... Args>
    friend SharedPtr<U> allocateShared(Allocator, Args&&...);

    template<typename U, typename... Args>
    friend SharedPtr<U> makeShared(Args&&...);

public:
    SharedPtr(): data(nullptr), cptr(nullptr) {}

    template<typename U, typename Deleter, typename Allocator>
    SharedPtr(U* data, Deleter deleter, Allocator allocator): data(data) {
        using ControlBlockType = ControlBlock<U, Deleter, Allocator>;
        using ControlBlockAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<ControlBlockType>;

        ControlBlockAllocator alloc(allocator);
        cptr = std::allocator_traits<ControlBlockAllocator>::allocate(alloc, 1);
        new(cptr) ControlBlockType(deleter, alloc);
    }

private:
    template<typename U, typename Deleter>
    void allocate_memory() {
        using ControlBlockType = ControlBlock<U, Deleter, std::allocator<U>>;

        std::allocator<ControlBlockType> alloc;
        cptr = std::allocator_traits<std::allocator<ControlBlockType>>::allocate(alloc, 1);
    }

public:
    template<typename U>
    explicit SharedPtr(U* data): data(data) {
        using ControlBlockType = ControlBlock<U, std::default_delete<U>, std::allocator<U>>;

        allocate_memory<U, std::default_delete<U>>();
        new(cptr) ControlBlockType();
    }

    template<typename U, typename Deleter>
    SharedPtr(U* data, Deleter deleter): data(data) {
        using ControlBlockType = ControlBlock<U, Deleter, std::allocator<U>>;

        allocate_memory<U, Deleter>();
        new(cptr) ControlBlockType(deleter);
    }

    template<typename U>
    SharedPtr(const SharedPtr<U>& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().shared_counter += 1;
        }
    }

    SharedPtr(const SharedPtr& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().shared_counter += 1;
        }
    }

    template <typename U>
    SharedPtr(SharedPtr<U>&& other):
            data(std::move(other.data)),
            cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    SharedPtr(SharedPtr&& other): data(std::move(other.data)), cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    template <typename U>
    SharedPtr(const WeakPtr<U>& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().shared_counter += 1;
        }
    }

    template <typename U, typename Allocator = std::allocator<U>>
    explicit SharedPtr(ControlBlock<U, std::default_delete<U>, Allocator>* cptr): cptr(cptr) {}

    T* get() const {
        if (data)  {
            return data;
        }
        if (cptr) {
            return &(static_cast<ControlBlock<T>*>(cptr)->data).value();
        }
        return nullptr;
    }

    void swap(SharedPtr& other) {
        std::swap(data, other.data);
        std::swap(cptr, other.cptr);
    }

    template<typename U>
    SharedPtr& operator=(const SharedPtr<U>& other) {
        SharedPtr save = other;
        swap(save);
        return *this;
    }

    SharedPtr& operator=(const SharedPtr& other) {
        SharedPtr save = other;
        swap(save);
        return *this;
    }

    template <typename U>
    SharedPtr& operator=(SharedPtr<U>&& other) {
        SharedPtr move_copy = std::move(other);
        swap(move_copy);
        return *this;
    }

    SharedPtr& operator=(SharedPtr&& other) {
        SharedPtr move_copy = std::move(other);
        swap(move_copy);
        return *this;
    }

    size_t use_count() const {
        if (cptr)
            return cptr->getCounter().shared_counter;
        return 0;
    }

    void reset() noexcept {
        SharedPtr().swap(*this);
    }

    template <typename U>
    void reset(U* data) {
        SharedPtr(data).swap(*this);
    }

    template <typename U, typename Deleter>
    void reset(U* ptr, Deleter deleter) {
        SharedPtr<T>(ptr, deleter).swap(*this);
    }

    template <typename U, typename Deleter, typename Allocator>
    void reset(U* ptr, Deleter deleter, Allocator allocator) {
        SharedPtr<T>(ptr, deleter, allocator).swap(*this);
    }

    T& operator*() const {
        if (data)
            return *data;
        return (static_cast<ControlBlock<T>*>(cptr)->data).value();
    }

    T* operator->() const {
        if (data)
            return data;
        return &(static_cast<ControlBlock<T>*>(cptr)->data).value();
    }

    ~SharedPtr() {
        if (cptr == nullptr) return;

        --(cptr->getCounter().shared_counter);
        if (use_count() > 0) return;

        if (data) {
            cptr->destroy_data(data);
        } else {
            cptr->destroy_partial();
        }

        if (cptr->getCounter().weak_counter == 0) {
            cptr->destroy_all();
        }
    }
};


template <typename T, typename Allocator = std::allocator<T>, typename... Args>
SharedPtr<T> allocateShared(Allocator alloc, Args&&... args) {
    using ControlBlockType = ControlBlock<T, std::default_delete<T>, Allocator>;
    using ControlBlockAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<ControlBlockType>;
    using AllocatorTraits = std::allocator_traits<ControlBlockAllocator>;

    ControlBlockAllocator cb_alloc;
    ControlBlockType* memory = AllocatorTraits::allocate(cb_alloc, 1);
    AllocatorTraits::construct(cb_alloc, memory, alloc, std::forward<Args>(args)...);
    memory->need_destuction = true;
    return SharedPtr<T>(memory);
}

template <typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
    return allocateShared<T>(std::allocator<T>(), std::forward<Args>(args)...);
}


template <typename T>
class WeakPtr {
private:
    T* data = nullptr;
    ControlBlockBase* cptr = nullptr;

public:
    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

    WeakPtr(): data(nullptr), cptr(nullptr) {}

    template <typename U>
    WeakPtr(const SharedPtr<U>& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().weak_counter += 1;
        }
    }

    WeakPtr(const SharedPtr<T>& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().weak_counter += 1;
        }
    }

    template <typename U>
    WeakPtr(SharedPtr<U>&& other):
            data(std::move(other.data)),
            cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    WeakPtr(SharedPtr<T>&& other):
            data(std::move(other.data)),
            cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    template <typename U>
    WeakPtr(const WeakPtr<U>& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().weak_counter += 1;
        }
    }

    WeakPtr(const WeakPtr& other): data(other.data), cptr(other.cptr) {
        if (cptr) {
            cptr->getCounter().weak_counter += 1;
        }
    }

    template <typename U>
    WeakPtr(WeakPtr<U>&& other): data(std::move(other.data)), cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    WeakPtr(WeakPtr&& other): data(std::move(other.data)), cptr(std::move(other.cptr)) {
        other.data = nullptr;
        other.cptr = nullptr;
    }

    void swap(WeakPtr& other) {
        std::swap(data, other.data);
        std::swap(cptr, other.cptr);
    }

    template<typename U>
    WeakPtr& operator=(const SharedPtr<U>& other) {
        WeakPtr save = other;
        swap(save);
        return *this;
    }

    WeakPtr& operator=(const SharedPtr<T>& other) {
        WeakPtr save = other;
        swap(save);
        return *this;
    }

    template<typename U>
    WeakPtr& operator=(SharedPtr<U>&& other) {
        WeakPtr save = std::move(other);
        swap(save);
        return *this;
    }

    WeakPtr& operator=(SharedPtr<T>&& other) {
        WeakPtr save = std::move(other);
        swap(save);
        return *this;
    }

    template<typename U>
    WeakPtr& operator=(const WeakPtr<U>& other) {
        WeakPtr save = other;
        swap(save);
        return *this;
    }

    WeakPtr& operator=(const WeakPtr<T>& other) {
        WeakPtr save = other;
        swap(save);
        return *this;
    }

    template<typename U>
    WeakPtr& operator=(WeakPtr<U>&& other) {
        WeakPtr save = std::move(other);
        swap(save);
        return *this;
    }

    WeakPtr& operator=(WeakPtr<T>&& other) {
        WeakPtr save = std::move(other);
        swap(save);
        return *this;
    }

    size_t use_count() const {
        if (cptr) {
            return cptr->getCounter().shared_counter;
        }
        return 0;
    }

    bool expired() const {
        return use_count() == 0;
    }

    SharedPtr<T> lock() const {
        if (use_count() == 0)
            return SharedPtr<T>();

        return SharedPtr<T>(*this);
    }

    ~WeakPtr() {
        if (cptr == nullptr) return;

        --(cptr->getCounter().weak_counter);
        if (cptr->getCounter().weak_counter == 0 && cptr->getCounter().shared_counter == 0) {
            cptr->destroy_all();
        }
    }
};