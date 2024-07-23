import io
import png
import mido
from mido import Message
from dataclasses import dataclass
from lib.pythonbmp.BITMAPlib import(
        newBMP, #insert other func here
        getfuncmetastr as meta,
        saveBMP,
        getcolorname2RGBdict,
        plotstring,
        copyrect,
        rectangle,
        convertselection2BMP,
        font8x14,
        font8x8,
        getBMPimgbytes
        )

@dataclass
class Chunk:
    x1: int
    y1: int
    x2: int
    y2: int
    def getSize(self)-> tuple[int:int]:
        return (self.x2 - self.x1, self.y2 - self.y1)

class ScreenData:
    def __init__(self, width, height, msg: mido.Message) -> None:
        self.width: int = width
        self.height: int = height
        m_bytes = msg.bytes()
        self.top_left: str = f"{hex(m_bytes[0] & 0xf0)}"
        self.top_right: str = f"{hex(m_bytes[0] & 0x0f)}"
        self.middle_left: str = f"{m_bytes[1]}"
        # self.value = m_bytes[2]
        self.middle_right: str = f"{m_bytes[2]}"
       

 
chunks: list[Chunk] = [
    Chunk(0  ,0  ,60, 60),
    Chunk(0  ,60 ,60, 80),
    Chunk(60 ,0  ,120, 60),
    Chunk(60 ,60  ,120, 80),
    Chunk(120 ,0 ,160, 60),
    Chunk(120 ,60 ,160, 80),
]

def least_sig_byte(value: int) -> tuple:
    # Extract the least significant two bytes
    ls_bytes = value & 0xFFFF
    # Split into two uint8 values
    byte2 = ls_bytes & 0xFF
    byte1 = (ls_bytes >> 8) & 0xFF
    return byte1, byte2

def magic_number(encoded_buffer_length: int) -> int:
    return (encoded_buffer_length // 128) * 128 - 8

def generate_sysex_header(chunk: Chunk, png_buffer_length: int, encoded_buffer_length: int) -> bytearray:
    sysex_array = bytearray([0x47, 0x7f, 0x4a, 0x04])

    message_size = encoded_buffer_length + 16
    high_byte, low_byte = least_sig_byte(message_size + magic_number(encoded_buffer_length))
    sysex_array.extend([high_byte, low_byte])

    low_len, high_len = least_sig_byte(png_buffer_length)
    if high_len >= 128:
        sysex_array.extend([0x20, 0x20])
        low_len, high_len = least_sig_byte(png_buffer_length - 128)
    else:
        sysex_array.extend([0x00, 0x20])

    sysex_array.extend([chunk.x1, 0x00, chunk.y1, 0x00])
    sysex_array.extend([high_len, low_len])

    return sysex_array


def encode_png(original_png_buffer: bytes) -> bytes:
    encode_png_buffer = bytearray()
    stride = 7

    i = 0
    while i < len(original_png_buffer):
        end = i + stride
        if end > len(original_png_buffer):
            end = len(original_png_buffer)
        group = bytearray(original_png_buffer[i:end])
        control_byte = 0

        for j in range(len(group)):
            if group[j] >= 128:
                control_byte |= 1 << j
                group[j] -= 128

        encode_png_buffer.append(control_byte)
        encode_png_buffer.extend(group)

        i += stride
    return bytes(encode_png_buffer)

def convert_bmp_to_png(width: int, height: int, bmp_bytes_array: list) -> bytes:
    image_array = [x for x in range(height)]
    row = []
    row_index = 0
    stride = (width * 3)
    i = 0
    for i in range(len(bmp_bytes_array)):
        color_index = 2 - (i % 3)

        pixel_count = i % stride
        row.append(bmp_bytes_array[i + color_index])
        if pixel_count == stride - 1:
            # image_array.append(row[:])
            image_array[(height-1) - row_index] = row[:]
            row.clear()
            row_index += 1

    image_buffer = io.BytesIO()
    w = png.Writer(width, height)
    w.write(image_buffer, image_array)
    return image_buffer.getvalue()

def generate_png_image_chunks(width: int, height: int, data: ScreenData) -> list[list]:
    bmp = newBMP(width, height, 24) # (x,y,bit depth)
    c = getcolorname2RGBdict() #friendly color names 2 rgb
    fontsize = 2 # font size
    pixspace = 0 # space between bitmap font pixels (0 = default)
    charspace = 0 # space bitmap font characters (0 = default)
    rectangle(bmp, 0, 0, width-1, height-1, c['black'])
    plotstring(bmp, 1, 5, data.top_left,fontsize, pixspace, charspace, c['brightyellow'], font8x8)
    plotstring(bmp, 90, 5, data.top_right,fontsize, pixspace, charspace, c['brightyellow'], font8x8)
    plotstring(bmp, 1, 25, data.middle_left,fontsize, pixspace, charspace, c['brightyellow'], font8x8)
    plotstring(bmp, 90, 25, data.middle_right,fontsize, pixspace, charspace, c['brightyellow'], font8x8)
    # rectangle(bmp, 1, 50, 80, 59, c['red'])

    image_chunks = []
    for chunk in chunks:
        buf = getBMPimgbytes(copyrect( bmp, chunk.x1, chunk.y1, chunk.x2 -1, chunk.y2-1 ))
        
        (cWidth, cHeight) = chunk.getSize()
        png_bytes = convert_bmp_to_png(cWidth, cHeight, buf[7:] )
        encoded_png_bytes = encode_png(png_bytes)
        image_chunks.append({"chunk": chunk, "original_png_bytes":png_bytes, "encoded_png_bytes": encoded_png_bytes})
    return image_chunks

def send_sysex_msg(out_port: any, chunk: Chunk, png_buffer_length: int, encoded_buffer_length: int, payload: bytes):
    sysex_array = generate_sysex_header(chunk, png_buffer_length, encoded_buffer_length)
    sysex_array.extend(payload)

    try:
        mMsg = mido.Message('sysex', data=sysex_array)
        out_port.send(mMsg)
    except Exception as e:
        print(f"ERROR: {e}")


def choose_output(midi_port_name: str):
    outPorts = mido.get_output_names()
    found_port = False
    outPortIndex: int
    for i, port in enumerate(outPorts):
        if port == midi_port_name:
            found_port = True
            outPortIndex = i
            break
    if not found_port:
        for i, port in enumerate(outPorts):
            print(f"{i}:, {port}")
        portNumber = input("Select output index ")
        outPortIndex = int(portNumber)
    return mido.open_output(outPorts[outPortIndex])

def choose_input(midi_port_name: str):
    in_ports = mido.get_input_names()
    found_port = False
    in_port_index: int
    for i, port in enumerate(in_ports):
        if port == midi_port_name:
            found_port = True
            in_port_index = i
            break
    if not found_port:
        for i, port in enumerate(in_ports):
            print(f"{i}:, {port}")
        in_port_selection = input("Select output index ")
        in_port_index = int(in_port_selection)
    return mido.open_input(in_ports[in_port_index])


def main():
    screen_width: int = 160
    screen_height: int = 80
    midi_port_name = "MPC Studio mk2 Public"
    in_port = choose_input(midi_port_name)
    out_port = choose_output(midi_port_name)
    msg: mido.Message
    for msg in in_port:
        screen_data = ScreenData(screen_width, screen_height, msg)
        image_chunks = generate_png_image_chunks(screen_width, screen_height, screen_data)
        for v in image_chunks:
            send_sysex_msg(out_port, v["chunk"], len(v["original_png_bytes"]), len(v["encoded_png_bytes"]), v["encoded_png_bytes"])
if __name__=="__main__":
        main()
