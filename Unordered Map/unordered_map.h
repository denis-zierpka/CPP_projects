#include <iostream>
#include <type_traits>
#include <vector>
#include <cmath>
#include <map>


template<typename T, typename Allocator = std::allocator<T>>
class List {
public:
    struct Node {
        T data;
        Node* next;
        Node* prev;

        Node() {}
        Node(const T& data): data(data) {}
        Node(const Node* obj): data(obj->data) {}
        Node(const Node& obj): data(obj.data) {}
        Node(T&& data): data(std::move(data)) {}
    };

    Node* head;
    size_t sz;

    using N_Allocator = typename Allocator::template rebind<Node>::other;
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
            return ptr != iter.ptr;
        }

        reference operator*() const {
            return ptr->data;
        }

        pointer operator->() const {
            return &(ptr->data);
        }
    };

    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    iterator insert_object(Node* object, const T& val) {
        Node* new_n = AllocTraits::allocate(n_allocator, 1);
        AllocTraits::construct(n_allocator, new_n, val);
        new_n->next = object->next;
        new_n->prev = object;
        object->next->prev = new_n;
        object->next = new_n;
        ++sz;
        iterator ans(new_n);
        return ans;
    }

    void insert_just_object(Node* object) {
        Node* new_n = AllocTraits::allocate(n_allocator, 1);
        AllocTraits::construct(n_allocator, new_n);
        new_n->next = object->next;
        new_n->prev = object;
        object->next->prev = new_n;
        object->next = new_n;
        ++sz;
    }

    iterator erase_object(Node* object) {
        iterator ans(object->next);
        object->prev->next = object->next;
        object->next->prev = object->prev;
        AllocTraits::destroy(n_allocator, object);
        AllocTraits::deallocate(n_allocator, object, 1);
        --sz;
        return ans;
    }

    explicit List(const Allocator& t_allocator = Allocator()): head(nullptr), sz(0), t_allocator(t_allocator) {
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
    }

    explicit List(size_t count, const T& data, const Allocator& t_allocator = Allocator()):
            head(nullptr),
            sz(0),
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
            sz(0),
            t_allocator(t_allocator) {
        head = AllocTraits::allocate(n_allocator, 1);
        head->next = head;
        head->prev = head;
        for (size_t i = 0; i < count; ++i) {
            insert_just_object(head->prev);
        }
    }

    List(const List& other): head(nullptr), sz(0) {
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

    List(List&& other): head(other.head), sz(other.sz) {
        auto new_h = std::allocator_traits<N_Allocator>::allocate(n_allocator, 1);
        new_h->next = new_h;
        new_h->prev = new_h;
        other.sz = 0;
        other.head = new_h;

        t_allocator = std::move(other.t_allocator);
        n_allocator = std::move(other.n_allocator);
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

        sz = 0;
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

    List& operator=(List&& other) noexcept {
        if (this == &other) return *this;

        List save = std::move(other);
        std::swap(sz, save.sz);
        std::swap(head, save.head);
        if (std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value && t_allocator != other.t_allocator) {
            t_allocator = std::move(other.t_allocator);
            n_allocator = std::move(other.n_allocator);
        }
        return *this;
    }


    size_t size() const {
        return sz;
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

    iterator insert(const_iterator ptr, const T& value) {
        return insert_object(const_cast<Node*>((--ptr).ptr), value);
    }

    iterator erase(const_iterator ptr) {
        return erase_object(const_cast<Node*>(ptr.ptr));
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

    void clear() {
        while (begin() != end()) {
            erase(begin());
        }
    }


    ~List() {
        Node* current = head->next;
        for (size_t i = 0; i < sz; ++i) {
            Node* save = current;
            current = current->next;
            AllocTraits::destroy(n_allocator, save);
            AllocTraits::deallocate(n_allocator, save, 1);
        }
        AllocTraits::deallocate(n_allocator, head, 1);
    }
};


template<
        typename Key,
        typename Value,
        typename Hash = std::hash<Key>,
        typename Equal = std::equal_to<Key>,
        typename Alloc = std::allocator<std::pair<const Key, Value>>
>
class UnorderedMap {
public:
    using NodeType = std::pair<const Key, Value>;

    using list_iter = typename List<NodeType*, typename Alloc::template rebind<NodeType*>::other>::iterator;
    using list_const_iter = typename List<NodeType*,typename Alloc::template rebind<NodeType*>::other>::const_iterator;
    template<bool IsConst>
    class CommonIterator {
    public:
        std::conditional_t<IsConst, list_const_iter, list_iter> ptr;

        using iterator_category = std::forward_iterator_tag;
        using value_type = typename std::conditional_t<IsConst, const NodeType, NodeType>;
        using reference = typename std::conditional_t<IsConst, const NodeType&, NodeType&>;
        using pointer = typename std::conditional_t<IsConst, const NodeType*, NodeType*>;
        using difference_type = std::ptrdiff_t;


        operator CommonIterator<true>() {
            return CommonIterator<true>(ptr);
        }

        CommonIterator(const std::conditional_t<IsConst, list_const_iter, list_iter>& node): ptr(node) {}

        CommonIterator& operator++() {
            ++ptr;
            return *this;
        }

        CommonIterator& operator--() {
            --ptr;
            return *this;
        }

        CommonIterator operator++(int) {
            auto save = *this;
            ++ptr;
            return save;
        }

        CommonIterator operator--(int) {
            auto save = *this;
            --ptr;
            return save;
        }

        bool operator==(const CommonIterator& other) const {
            return ptr == other.ptr;
        }

        bool operator!=(const CommonIterator& other) const {
            return ptr != other.ptr;
        }

        reference operator*() const {
            return **ptr;
        }

        pointer operator->() const {
            return *ptr;
        }
    };


    using Iterator = CommonIterator<false>;
    using ConstIterator = CommonIterator<true>;



    using NAlloc = typename Alloc::template rebind<NodeType>::other;

    using AllocTraits = std::allocator_traits<NAlloc>;
    NAlloc n_allocator;
    List<NodeType*, typename Alloc::template rebind<NodeType*>::other> elements;
    Hash hashfn;
    Equal equalfn;
    std::vector<list_iter> hash_table;
    std::vector<list_iter> save_last_hash;

    static const int resize_value = 1;
    float max_load_factor_value = 0.75;

    UnorderedMap() {
        hash_table.resize(resize_value, elements.end());
        save_last_hash.resize(resize_value, elements.end());
    }

    explicit UnorderedMap(const Alloc& alloc): n_allocator(alloc) {
        hash_table.resize(resize_value, elements.end());
        save_last_hash.resize(resize_value, elements.end());
    }

    void move_from_other(UnorderedMap&& other) {
        elements = std::move(other.elements);
        hashfn = std::move(other.hashfn);
        equalfn = std::move(other.equalfn);
        hash_table = std::move(other.hash_table);
        save_last_hash = std::move(other.save_last_hash);
    }

    void get_from_other(const UnorderedMap& other) {
        hashfn = other.hashfn;
        equalfn = other.equalfn;
    }

    UnorderedMap(const UnorderedMap& other) {
        n_allocator = AllocTraits::select_on_container_copy_construction(other.n_allocator);
        get_from_other(other);
        hash_table.resize(resize_value, elements.end());
        save_last_hash.resize(resize_value, elements.end());
        for (auto& ptr : other) {
            insert(ptr);
        }
    }

    UnorderedMap(UnorderedMap&& other) {
        n_allocator = std::move(AllocTraits::select_on_container_copy_construction(other.n_allocator));
        move_from_other(std::move(other));
    }

    UnorderedMap& operator=(const UnorderedMap& other) {
        if (this == &other) return *this;

        n_allocator = other.n_allocator;
        get_from_other(other);

        hash_table.resize(resize_value, elements.end());
        save_last_hash.resize(resize_value, elements.end());
        for (auto& ptr : other) {
            insert(ptr);
        }
        return *this;
    }

    UnorderedMap& operator=(UnorderedMap&& other) noexcept {
        if (this == &other) return *this;

        n_allocator = std::move(other.n_allocator);
        move_from_other(std::move(other));
        return *this;
    }

    size_t bucket_pos(const Key& key) const {
        return hashfn(key) % hash_table.size();
    }

    Alloc get_allocator() const noexcept {
        return n_allocator;
    }

    size_t size() const {
        return elements.size();
    }

    Iterator begin() {
        return static_cast<Iterator>(elements.begin());
    }

    ConstIterator begin() const {
        return static_cast<ConstIterator>(elements.cbegin());
    }

    ConstIterator cbegin() const {
        return static_cast<ConstIterator>(elements.cbegin());
    }

    Iterator end() {
        return static_cast<Iterator>(elements.end());
    }

    ConstIterator end() const {
        return static_cast<ConstIterator>(elements.cend());
    }

    ConstIterator cend() const {
        return static_cast<ConstIterator>(elements.cend());
    }

    size_t max_size() const noexcept {
        return std::numeric_limits<size_t>::max();
    }

    Iterator find(const Key& key, size_t pos = 0, bool pos_calculated = false) {
        if (!pos_calculated)
            pos = bucket_pos(key);
        list_iter ptr = hash_table[pos];

        while (ptr != elements.end()) {
            if (equalfn((*ptr)->first, key)) {
                return static_cast<Iterator>(ptr);
            }
            if (ptr == save_last_hash[pos]) {
                break;
            }
            ++ptr;
        }
        return static_cast<Iterator>(end());
    }

    float load_factor() const {
        return static_cast<float>(elements.size()) / hash_table.size();
    }

    void max_load_factor(float value) {
        max_load_factor_value = value;
    }

    float max_load_factor() {
        return max_load_factor_value;
    }

    void check_load_factor() {
        if (load_factor() > max_load_factor()) {
            reserve(hash_table.size() * 2 + 1);
        }
    }

    void reserve(size_t count) {
        if (count > hash_table.size()) {
            rehash(std::ceil(count / max_load_factor()));
        }
    }

    void rehash(size_t count) {
        hash_table.clear();
        save_last_hash.clear();
        auto save = std::move(elements);
        hash_table.resize(count, elements.end());
        save_last_hash.resize(count, elements.end());
        for (list_iter i = save.begin(); i != save.end(); ++i) {
            list_iter& elem = hash_table[bucket_pos((*i)->first)];
            bool update_last = false;
            if (elem == elements.end()) {
                update_last = true;
            }

            elem = elements.insert(elem, *i);
            if (update_last) {
                save_last_hash[bucket_pos((*i)->first)] = elem;
            }
        }
    }


    template<typename ...Args>
    std::pair<Iterator, bool> emplace(Args&&... args) {
        NodeType* ptr = AllocTraits::allocate(n_allocator, 1);
        AllocTraits::construct(n_allocator, ptr, std::forward<Args>(args)...);

        Iterator res = find(ptr->first);
        if (res != end())
            return {res, false};

        check_load_factor();
        size_t save_bucket_pos = bucket_pos(ptr->first);
        list_iter save = hash_table[save_bucket_pos];
        list_iter ans = elements.insert(save, ptr);

        if (elements.end() == save) {
            hash_table[save_bucket_pos] = ans;
            save_last_hash[save_bucket_pos] = ans;

        }
        return {static_cast<Iterator>(ans), true};
    }

    std::pair<Iterator, bool> insert(const NodeType& value) {
        return emplace(value);
    }

    std::pair<Iterator, bool> insert(NodeType&& value) {
        return emplace(std::move(value));
    }

    template<class P>
    std::pair<Iterator, bool> insert(P&& value) {
        return emplace(std::forward<P>(value));
    }

    template<typename InputIt>
    void insert(InputIt first, InputIt last) {
        for (InputIt i = first; i != last; ++i) {
            insert(*i);
        }
    }


    size_t erase(const Key& key) {
        std::pair<Iterator, size_t> ptr = r_find(key);
        if (ptr.first == elements.end()) {
            return 0;
        }
        r_erase(ptr.first, ptr.second, true);
        return 1;
    }


    Iterator erase(ConstIterator ptr) {
        return r_erase(ptr);
    }

    Iterator erase(ConstIterator first, ConstIterator last) {
        while (first != last) {
            first = erase(first);
        }
        using ListSave = typename List<NodeType*, typename Alloc::template rebind<NodeType*>::other>::Node*;
        return Iterator(const_cast<ListSave>(first.ptr.ptr));
    }

    Value& at(const Key& key) {
        Iterator ptr = find(key);
        if (ptr != end()) {
            return ptr->second;
        }
        throw std::out_of_range("Out of range!");
    }

    Value& operator[](const Key& key) {
        Iterator ptr_check = find(key);
        if (ptr_check != end()) {
            return at(key);
        }
        Iterator ptr = insert({key, Value()}).first;
        return ptr->second;
    }

    ~UnorderedMap() {
        for (list_const_iter ptr = elements.begin(); ptr != elements.end(); ++ptr) {
            delete *ptr;
        }
    }

private:
    std::pair<Iterator, size_t> r_find(const Key& key) {
        size_t save_bucket_pos = bucket_pos(key);
        return {find(key, save_bucket_pos, true), save_bucket_pos};
    }

    Iterator r_erase(ConstIterator ptr, size_t hash = 0, bool hash_calculated = false) {
        if (!hash_calculated)
            hash = bucket_pos(ptr->first);
        auto del = elements.erase(ptr.ptr);
        if (static_cast<ConstIterator>(hash_table[hash]) == ptr) {
            if (del != elements.end() && bucket_pos((*del)->first) == hash)
                hash_table[hash] = del;
            else
                hash_table[hash] = elements.end();
            return static_cast<Iterator>(del);
        }
        return static_cast<Iterator>(del);
    }

};

