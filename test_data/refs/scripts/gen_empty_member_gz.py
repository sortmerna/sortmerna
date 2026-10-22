"""Generate multi-member FASTQ.gz with empty members inserted (issue #466 repro).
usage: gen.py in.fastq out.gz copies reads_per_member empty_kind(stored|fixed|none) [every]"""
import sys, zlib, struct
src, out, copies, per, kind = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
every = int(sys.argv[6]) if len(sys.argv) > 6 else 37
lines = open(src, 'rb').read().splitlines(keepends=True)
recs = [b''.join(lines[i:i+4]) for i in range(0, len(lines), 4)]
HDR = bytes.fromhex('1f8b08000000000000ff')
EMPTY = {'stored': HDR + bytes.fromhex('010000ffff') + b'\0'*8,
         'fixed':  HDR + bytes.fromhex('0300') + b'\0'*8, 'none': b''}[kind]
def member(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    return HDR + c.compress(data) + c.flush() + struct.pack('<II', zlib.crc32(data), len(data) & 0xffffffff)
n = 0; nempty = 0
with open(out, 'wb') as f:
    buf = []
    for k in range(copies):
        for r in recs:
            buf.append(r); n += 1
            if len(buf) == per:
                f.write(member(b''.join(buf))); buf = []
                if (n // per) % every == 0:   # an empty member every `every` members
                    f.write(EMPTY); nempty += 1
    if buf: f.write(member(b''.join(buf)))
print(out, 'reads', n, 'empty members', nempty if kind != 'none' else 0)
