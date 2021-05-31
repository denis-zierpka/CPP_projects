#include <iostream>
#include <vector>
#include <type_traits>


template<size_t chunkSize>
class FixedAllocator {
private:
    std::vector<char*> free_positions_;
    static const size_t kFreeCount = 40;

    static FixedAllocator<chunkSize>& instance() noexcept {
        static auto item = FixedAllocator<chunkSize>();
        return item;
    };

    template<typename>
    friend class FastAllocator;

public:
    void* allocate() {
        if (free_positions_.empty()) {
            char* new_pool = reinterpret_cast<char*>(::operator new(chunkSize * kFreeCount));
            for (size_t i = 0; i < kFreeCount; ++i) {
                free_positions_.push_back(new_pool + i * chunkSize);
            }
        }
        char* result = free_positions_.back();
        free_positions_.pop_back();
        return reinterpret_cast<void*>(result);
    }

    void deallocate(void* ptr) {
        free_positions_.push_back(reinterpret_cast<char*>(ptr));
    }
};


template<typename T>
class FastAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;

    template<typename U>
    explicit FastAllocator(const FastAllocator<U>&) {}

    FastAllocator() = default;

    T* allocate(size_t count) {
        size_t result_count = count * sizeof(T);

        if (result_count <= 4) {
            return reinterpret_cast<T*>(FixedAllocator<4>::instance().allocate());
        }
        if (result_count <= 8) {
            return reinterpret_cast<T*>(FixedAllocator<8>::instance().allocate());
        }
        if (result_count <= 16) {
            return reinterpret_cast<T*>(FixedAllocator<16>::instance().allocate());
        }
        if (result_count <= 24) {
            return reinterpret_cast<T*>(FixedAllocator<24>::instance().allocate());
        }
        return reinterpret_cast<T*>(::operator new(count * sizeof(T)));
    }

    void deallocate(T* ptr, size_t count) {
        size_t result_count = count * sizeof(T);

        if (result_count <= 4) {
            FixedAllocator<4>::instance().deallocate(reinterpret_cast<void*>(ptr));
            return;
        }
        if (result_count <= 8) {
            FixedAllocator<8>::instance().deallocate(reinterpret_cast<void*>(ptr));
            return;
        }
        if (result_count <= 16) {
            FixedAllocator<16>::instance().deallocate(reinterpret_cast<void*>(ptr));
            return;
        }
        if (result_count <= 24) {
            FixedAllocator<24>::instance().deallocate(reinterpret_cast<void*>(ptr));
            return;
        }
        ::operator delete(ptr);
    }
};


template<typename T, typename U>
bool operator==(const FastAllocator<T>& a, const FastAllocator<U>& b) {
    return true;
}

template<typename T, typename U>
bool operator!=(const FastAllocator<T>& a, const FastAllocator<U>& b) {
    return !(::operator==(a, b));
}


template<typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct Node {
        T data;
        Node* next;
        Node* prev;

        Node() {}
        Node(const T& data): data(data) {}
        Node(const Node* obj): data(obj->data) {}
        Node(const Node& obj): data(obj.data) {}
    };
    Node* head;
    size_t size_;
    
    using Alloc_AllocTraits = std::allocator_traits<Allocator>;
    using N_Allocator = typename Alloc_AllocTraits::template rebind_alloc<Node>;
    using AllocTraits = std::allocator_traits<N_Allocator>;
    N_Allocator n_allocator;
    Allocator t_allocator;

    template<bool IsConst>
    class Iterator {
    public:
        std::conditional_t<IsConst, const Node*, Node*> ptr;

        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = typename std::conditional<IsConst, const T, T>::type;
        using reference = typename std::conditional<IsConst, const T&, T&>::type;
        using pointer = typename std::conditional<IsConst, const T*, T*>::type;
        using difference_type = std::ptrdiff_t;


        operator Iterator<true>() {
            return Iterator<true>(ptr);
        }

        Iterator(const Iterator& iter): ptr(iter.ptr) {}

        Iterator(std::conditional_t<IsConst, const Node*, Node*> node_ptr): ptr(node_ptr) {}

        Iterator& operator++() {
            ptr = ptr->next;
            return *this;
        }

        Iterator& operator--() {
            ptr = ptr->prev;
            return *this;
        }

        Iterator operator++(int) {
            Iterator temp = *this;
            ptr = ptr->next;
            return temp;
        }

        Iterator operator--(int) {
            Iterator temp = *this;
            ptr = ptr->prev;
            return temp;
        }

        bool operator==(const Iterator& iter) const {
            return ptr == iter.ptr;
        }

        bool operator!=(const Iterator& iter) const {
            return !(*this == iter);
        }

        reference operator*() const {
            return ptr->data;
        }

        pointer operator->() const {
            return &(ptr->data);
        }
    };

public:
    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    void link_to_object(Node* first, Node* second) {
        second->next = first->next;
        second->prev = first;
        first->next->prev = second;
        first->next = second;
    }

    void insert_object(Node* object, const T& val) {
        Node* new_n = AllocTraits::allocate(n_allocator, 1);
        AllocTraits::construct(n_allocator, new_n, val);
        link_to_object(object, new_n);
        ++size_;
    }

    void insert_just_object(Node* object) {
        Node* new_n = AllocTraits::allocate(n_allocator, 1);
        AllocTraits::construct(n_allocator, new_n);
        link_to_object(object, new_n);
        ++size_;
    }

    void erase_object(Node* object) {
        object->prev->next = object->next;
        object->next->prev = object->prev;
        AllocTraits::destroy(n_allocator, object);
        AllocTraits::deallocate(n_allocator, object, 1);
        --size_;
    }

public:
    explicit List(const Allocator& t_allocator = Allocator()): head(nullptr), size_(0), t_allocator(t_allocator) {
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
    }

    explicit List(size_t count, const T& data, const Allocator& t_allocator = Allocator()):
            head(nullptr),
            size_(0),
            t_allocator(t_allocator) {
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
        for (size_t i = 0; i < count; ++i) {
            push_back(data);
        }
    }

    explicit List(size_t count, const Allocator& t_allocator = Allocator()):
            head(nullptr),
            size_(0),
            t_allocator(t_allocator) {
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
        for (size_t i = 0; i < count; ++i) {
            insert_just_object(head->prev);
        }
    }

    List(const List& other): head(nullptr), size_(0) {
        t_allocator = std::allocator_traits<Allocator>::select_on_container_copy_construction(other.t_allocator);
        n_allocator = std::allocator_traits<N_Allocator>::select_on_container_copy_construction(other.n_allocator);
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
        Node* other_ptr = other.head;
        for (size_t i = 0; i < other.size(); ++i) {
            push_back(other_ptr->next->data);
            other_ptr = other_ptr->next;
        }
    }

    List& operator=(const List& other) {
        if (this == &other) return *this;

        Node* current = head->next;
        for (size_t i = 0; i < size(); ++i) {
            Node* save = current;
            current = current->next;
            AllocTraits::destroy(n_allocator, save);
            AllocTraits::deallocate(n_allocator, save, 1);
        }
        AllocTraits::deallocate(n_allocator, head, 1);

        size_ = 0;
        head = nullptr;
        if (AllocTraits::propagate_on_container_copy_assignment::value) {
            t_allocator = other.t_allocator;
            n_allocator = other.n_allocator;
        }
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
        Node* other_ptr = other.head;
        for (size_t i = 0; i < other.size(); ++i) {
            push_back(other_ptr->next->data);
            other_ptr = other_ptr->next;
        }
        return *this;
    }

    size_t size() const {
        return size_;
    }

    Allocator get_allocator() const {
        return t_allocator;
    }

    void push_front(const T& value) {
        insert_object(head, value);
    }

    void push_back(const T& value) {
        insert_object(head->prev, value);
    }

    void pop_front() {
        erase_object(head->next);
    }

    void pop_back() {
        erase_object(head->prev);
    }

    void insert(const_iterator ptr, const T& value) {
        insert_object(const_cast<Node*>((--ptr).ptr), value);
    }

    void erase(const_iterator ptr) {
        erase_object(const_cast<Node*>(ptr.ptr));
    }

    iterator begin() const {
        iterator ptr(head->next);
        return ptr;
    }

    iterator end() const {
        iterator ptr(head);
        return ptr;
    }

    const_iterator cbegin() const {
        const_iterator ptr(head->next);
        return ptr;
    }

    const_iterator cend() const {
        const_iterator ptr(head);
        return ptr;
    }

    reverse_iterator rbegin() const {
        reverse_iterator ptr(head);
        return ptr;
    }

    reverse_iterator rend() const {
        reverse_iterator ptr(head->next);
        return ptr;
    }

    const_reverse_iterator crbegin() const {
        const_reverse_iterator ptr(head);
        return ptr;
    }

    const_reverse_iterator crend() const {
        const_reverse_iterator ptr(head->next);
        return ptr;
    }

    ~List() {
        Node* current = head->next;
        for (size_t i = 0; i < size_; ++i) {
            Node* save = current;
            current = current->next;
            AllocTraits::destroy(n_allocator, save);
            AllocTraits::deallocate(n_allocator, save, 1);
        }
        AllocTraits::deallocate(n_allocator, head, 1);
    }
};
