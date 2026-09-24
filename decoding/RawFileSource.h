#include "RawSource.h"

#include <fstream>

struct EventHeader {

    uint64_t event_id;
    uint64_t timestamp;
};

class RawFileSource : public RawSource {
public:
    ~RawFileSource() override;

    bool open(const std::string& path) override;
    std::optional<RawEvent> next() override;
    uint64_t events_read() const noexcept override { return events_read; }
    void close() override;

private:
    bool read_block(std::span<std::byte> block);

    bool read_bank(RawEvent& event);

    std::ifstream file;
    uint64_t      events_read = 0;

};
