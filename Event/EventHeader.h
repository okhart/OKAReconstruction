class EventHeader {
public:
    EventHeader() = default;
    ~EventHeader() = default;

    void setRunNumber(uint32_t number);
    void setEventNumber(uint64_t number);
    void setTimestamp(uint64_t timestamp);

    uint32_t getRunNumber() const;
    uint64_t getEventNumber() const;
    uint64_t getTimestamp() const;

private:
    uint32_t runNumber = 0;
    uint64_t eventNumber = 0;
    uint64_t timestamp = 0;
};