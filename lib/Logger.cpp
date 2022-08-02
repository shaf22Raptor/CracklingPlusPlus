// Logger Class and Log Buffer Class
#include "../include/Logger.hpp"

using std::string;
using std::ostream;
using std::ofstream;

logBuffer::logBuffer() = default;

logBuffer::logBuffer(const ostream& stdOut, const ostream& destFile) : fileBuf(destFile.rdbuf()), stdBuf(stdOut.rdbuf()) {}

int logBuffer::sync() {
    // Call sync for both output streams and buffers
    return stdBuf->pubsync() & fileBuf->pubsync();
}

int logBuffer::underflow(int c)
{
    // Do nothing, return EOF
    return EOF;
}

int logBuffer::overflow(int c)
{
    // Insert Timestamp
    if (isAtStartOfLine) {
        time_t rawtime = time(nullptr);
        struct tm timeinfo;
        p_localtime(timeinfo, rawtime);
        string timestampBuffer;
        timestampBuffer.reserve(32);
        strftime(&timestampBuffer[0], timestampBuffer.size(), ">>> %Y-%m-%d %H:%M:%S: ", &timeinfo);
        stdBuf->sputn(&timestampBuffer[0], timestampBuffer.size());
        fileBuf->sputn(&timestampBuffer[0], timestampBuffer.size());
    }
    // Check for new line
    isAtStartOfLine = c == '\n';
    // Add extra newline
    if (isAtStartOfLine) { stdBuf->sputc((char)c); fileBuf->sputc((char)c); }
    // Pass 'c' to output buffers
    return stdBuf->sputc((char)c) & fileBuf->sputc((char)c);
}

Logger::Logger(ostream& logSource, const string& outFile) : origOutputStream(&logSource), origSrcBuffer(logSource.rdbuf())
{
    // Open file stream
    outputFileStream = ofstream(outFile, std::ios::binary);

    // Create custom logBuffer
    customBuffer = logBuffer(*origOutputStream, outputFileStream);

    // Create new ostream
    ostream loggedOutputStream(&customBuffer);

    // Set log source to use new ostream that uses our custom buffer
    origOutputStream->rdbuf(loggedOutputStream.rdbuf());
}

void Logger::close() {
    // Close file stream
    outputFileStream.close();
    // Reset original ref
    origOutputStream->rdbuf(origSrcBuffer);
}