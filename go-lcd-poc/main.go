package main

import (
	"bufio"
	"bytes"
	"fmt"
	"image"
	"image/png"
	"os"

	"github.com/fogleman/gg"
	"gitlab.com/gomidi/midi/v2"
	"gitlab.com/gomidi/midi/v2/drivers"
	_ "gitlab.com/gomidi/midi/v2/drivers/rtmididrv"
	"gopkg.in/yaml.v3"
)

type MPCControl struct {
	Name    string `yaml:"name"`
	Status  byte   `yaml:"status"`
	Channel byte   `yaml:"channel"`
	Data1   byte   `yaml:"data1"`
}

type MPC struct {
	Controls []MPCControl `yaml:"controls"`
}

func (c *MPC) GetControlByMsg(msg *midi.Message) *MPCControl {
	bytes := msg.Bytes()
	status := bytes[0] & 0xf0
	statuses := []byte{0x80, 0x90, 0xA0}
	if status == statuses[0] || status == statuses[1] || status == statuses[2] {
		status = 0x90
	}
	channel := bytes[0] & 0x0f
	data1 := bytes[1]
	for i := range c.Controls {
		ctrl := c.Controls[i]
		if ctrl.Status == status && ctrl.Data1 == data1 && ctrl.Channel == channel {
			return &ctrl
		}
	}
	return nil
}

var screen_width int = 160
var screen_height int = 80
var chuckDef = [6]image.Rectangle{
	image.Rect(0, 0, 60, 60),
	image.Rect(0, 60, 60, 80),
	image.Rect(60, 0, 120, 60),
	image.Rect(60, 60, 120, 80),
	image.Rect(120, 0, 160, 60),
	image.Rect(120, 60, 160, 80),
}
var bufIndex = [6]*bytes.Buffer{
	new(bytes.Buffer),
	new(bytes.Buffer),
	new(bytes.Buffer),
	new(bytes.Buffer),
	new(bytes.Buffer),
	new(bytes.Buffer),
}

func getPorts(inPortName string, outPortName string) (drivers.In, drivers.Out, error) {
	in_port, err := midi.FindInPort(inPortName)
	if err != nil {
		return nil, nil, err
	}
	out_port, err := midi.FindOutPort(outPortName)
	if err != nil {
		return nil, nil, err
	}
	return in_port, out_port, nil
}

// Convert image.Image to []byte
func image_to_png_byte(chunk_index int, img image.Image) ([]byte, error) {
	buf := bufIndex[chunk_index]
	buf.Reset()
	encoder := png.Encoder{CompressionLevel: png.BestCompression}
	err := encoder.Encode(buf, img)
	if err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}
func draw_image(ctx *gg.Context, mpc *MPC, msg *midi.Message) {
	ctx.SetRGBA(0, 0, 0, 1)
	ctx.Clear()
	status := msg.Bytes()[0] & 0xf0
	data1 := msg.Bytes()[1]
	value := msg.Bytes()[2]
	start := 1
	end := int((float64(value) / float64(127)) * float64((screen_width - 5)))
	ctx.SetRGBA(0, 1, 0, 1)
	ctx.DrawRectangle(float64(start), 10, float64(end), 20)
	ctx.Fill()
	ctx.SetRGBA(0, 1, 0, 1)
	if err := ctx.LoadFontFace("OpenSans.ttf", 16); err != nil {
		panic(err)
	}
	ctrl := mpc.GetControlByMsg(msg)
	if ctrl != nil {

		ctx.DrawString(fmt.Sprintf("%v", ctrl.Name), float64(start), 50)
	}
	ctx.DrawString(fmt.Sprintf("[ %X, %X, %X ]", status, data1, value), float64(start), 70)

}

func magicNumber(encodedBufferLength int) int {
	return int(encodedBufferLength)/128*128 - 8
}
func leastSigByte(value int) (byte, byte) {
	// Extract the least significant two bytes
	lsBytes := value & 0xFFFF
	// Split into two uint8 values
	byte2 := uint8(lsBytes & 0xFF)
	byte1 := uint8((lsBytes >> 8) & 0xFF)
	return byte(byte1), byte(byte2)
}
func generate_sysex_header(chunk image.Rectangle, pngBufferLength int, encodedBufferLength int) []byte {
	// 1: F0 (1 byte) SysEx Message Start
	// 2: 47 7F 4A (3 bytes)	Manufacturer / Device / Model ID
	// 3: 04 (1 byte) Message ID
	sysex_array := []byte{0x47, 0x7f, 0x4a, 0x04}

	// 4: XX XX (2 bytes) Size of the message + a magic number,
	message_size := encodedBufferLength + 16
	high_byte, low_byte := leastSigByte(message_size + magicNumber(encodedBufferLength))
	sysex_array = append(sysex_array, high_byte, low_byte)

	// 5: 00 / 20 (1 byte)	If 20, add 128 to payload size (8)
	var high_len, low_len byte
	low_len, high_len = leastSigByte(pngBufferLength)
	if high_len >= 128 {
		sysex_array = append(sysex_array, byte(0x20), byte(0x20))
		low_len, high_len = leastSigByte(pngBufferLength - 128)
	} else {
		sysex_array = append(sysex_array, byte(0x00), byte(0x20))

	}
	// 7: XX 00 YY 00 (4 bytes)	X / Y coordinates of the chunk
	sysex_array = append(sysex_array, byte(chunk.Min.X), 0x00, byte(chunk.Min.Y), 0x00)

	// 8: XX XX (2 bytes)	Size of the image (decoded to png!), but the bytes are swapped (so 1 byte is encoded as 01 00
	sysex_array = append(sysex_array, high_len, low_len)
	return sysex_array

}

// send_sysex_msg: Sends a write screen sysex message to the MPC Studio Black.
func send_sysex_msg(out_port *drivers.Out, chunk image.Rectangle, pngBufferLength int, encodedBufferLength int, payload *[]byte) {
	sysex_array := generate_sysex_header(chunk, pngBufferLength, encodedBufferLength)
	sysex_array = append(sysex_array, *payload...)
	send, err := midi.SendTo(*out_port)
	if err != nil {
		fmt.Printf("ERROR: %s\n", err)
		return
	}
	mMsg := midi.SysEx(sysex_array)
	err = send(mMsg)
	if err != nil {
		fmt.Printf("ERROR: %s\n", err)
	}

}

func write_to_screen(out_port *drivers.Out, image_data *image.RGBA) error {

	for chunkIndex := range len(chuckDef) {

		chunk := chuckDef[chunkIndex]
		chunk_img := (*image_data).SubImage(chunk)
		chunk_buffer_bytes, err := image_to_png_byte(chunkIndex, chunk_img)
		if err != nil {
			fmt.Printf("ERROR: %s", err)
			return err
		}
		encodedBuffer := encodePNG(chunk_buffer_bytes)
		send_sysex_msg(out_port, chunk, len(chunk_buffer_bytes), len(encodedBuffer), &encodedBuffer)
	}
	return nil
}
func encodePNG(originalBuffer []byte) []byte {
	var encodedBuffer []byte
	stride := 7
	for i := 0; i < len(originalBuffer); i += stride {
		end := i + stride
		if end > len(originalBuffer) {
			end = len(originalBuffer)
		}
		group := originalBuffer[i:end]
		controlByte := 0

		for j := 0; j < len(group); j++ {
			if group[j] >= 128 {
				controlByte |= 1 << j
				group[j] -= 128
			}
		}
		encodedBuffer = append(encodedBuffer, byte(controlByte))
		encodedBuffer = append(encodedBuffer, group...)
	}
	return encodedBuffer
}
func list_midi_ports() {
	for inPortIndex := range midi.GetInPorts() {

		inPort, err := midi.InPort(inPortIndex)
		if err != nil {
			fmt.Printf("ERROR: %s\n", err)
		}
		fmt.Printf("%d In: %s\n", inPortIndex, inPort.String())
	}
	for outPortIndex := range midi.GetOutPorts() {

		inPort, err := midi.InPort(outPortIndex)
		if err != nil {
			fmt.Printf("ERROR: %s\n", err)
		}
		fmt.Printf("%d Out: %s\n", outPortIndex, inPort.String())
	}
}

func main() {
	dummyImage := image.NewRGBA(image.Rect(0, 0, screen_width, screen_height))
	newCtx := gg.NewContextForRGBA(dummyImage)
	var controls MPC
	dat, err := os.ReadFile("controls.yaml")
	if err != nil {
		fmt.Println(err.Error())
	}
	err = yaml.Unmarshal(dat, &controls)
	if err != nil {
		fmt.Printf("ERROR: %v", err.Error())
	}

	// List Midi Ports, In/Out
	// list_midi_ports()
	defer midi.CloseDriver()
	inPort, outPort, err := getPorts("MPC Studio mk2 Public", "MPC Studio mk2 Public")
	if err != nil {
		fmt.Println(err.Error())
	}
	fmt.Printf("Using Input: %s\n", inPort.String())
	fmt.Printf("Using Output: %s\n", outPort.String())
	stop, err := midi.ListenTo(inPort, func(msg midi.Message, timestampms int32) {

		// Generate image: Prints a blank screen with the string of midi msg. Outputs Image data in RGBA in uint32 (r,g,b,a) format
		draw_image(newCtx, &controls, &msg)
		err := write_to_screen(&outPort, dummyImage)
		if err != nil {
			fmt.Printf("ERROR: %s\n", err)
			return
		}

	}, midi.UseSysEx())

	if err != nil {
		fmt.Printf("ERROR: %s\n", err)
		return
	}

	reader := bufio.NewReader(os.Stdin)
	input, err := reader.ReadString('\n')
	if err != nil {
		fmt.Println(err.Error())
		return
	}
	fmt.Println(input)

	stop()
}
