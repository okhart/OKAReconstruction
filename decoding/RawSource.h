class RawSource {
public:
    virtual ~RawSource() = default;

    virtual bool open(const std::string& location) = 0;

    virtual std::optional<RawEvent> next() = 0;

    virtual uint64_t events_read() const noexcept = 0;

    virtual void close() = 0;
};
