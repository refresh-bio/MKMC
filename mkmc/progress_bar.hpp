#ifndef _PROGRESS_BAR_
#define _PROGRESS_BAR_

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/ioctl.h>
#endif

#include <iostream>
#include <string>
#include <mutex>
#include <atomic>


class ProgressBar {
  public:
    ProgressBar(uint64_t total,
                const std::string &description = "",
                std::ostream &out = std::cerr,
                bool silent = false);

    ~ProgressBar();

    uint64_t GetTotal() const { return total_; }
    void SetFrequencyUpdate(uint64_t frequency_update_);
    void SetStyle(char unit_bar, char unit_space);

    ProgressBar& operator++();
    ProgressBar& operator+=(uint64_t delta);

  private:
    ProgressBar(const ProgressBar &) = delete;
    ProgressBar& operator=(const ProgressBar &) = delete;

    void ShowProgress(uint64_t progress) const;
    int GetConsoleWidth() const;
    int GetBarLength() const;

    bool silent_;
    bool logging_mode_;
    uint64_t total_;
    std::atomic<uint64_t> progress_ = {0};
    uint64_t frequency_update;
    std::ostream *out;
    mutable std::mutex mu_;
    mutable std::string buffer_;

    std::string description_;
    char unit_bar_ = '=';
    char unit_space_ = ' ';
};

//++ on atomic for each element may be very costly
//this class may be used to limit number of this increments
//its not thread safe!
//the object must live longer than passed progress_bar
//because it uses it in dtor
class ProgressBarUpdater {
    ProgressBar& progress_bar;
    uint64_t frequency_update;
    uint64_t cur_val{};
public:
    ProgressBarUpdater(ProgressBar& progress_bar, uint64_t frequency_update):
        progress_bar(progress_bar), frequency_update(frequency_update) {

    }
    ProgressBarUpdater& operator++() {
        ++cur_val;
        if (cur_val == frequency_update) {
            progress_bar += cur_val;
            cur_val = 0;
        }
        return *this;
    }
    ~ProgressBarUpdater() {
        progress_bar += cur_val;
        cur_val = 0; //not needed
    }
};

#endif // _PROGRESS_BAR_
