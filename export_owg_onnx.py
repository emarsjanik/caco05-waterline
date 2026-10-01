#!/usr/bin/env python3
"""
Export A Trained Wave Gauge So The Station Can Run It Without TensorFlow
===========================================================================
Run in Colab, next to the weights. Writes <stem>.onnx: the trained network
as one portable file that OpenCV (already on the station) or onnxruntime can
run. The station has no TensorFlow, and its Python 3.8 is too old for the
TensorFlow the model was trained with, so the weights file alone is no use
there.

CHECK FOR EQUAL ANSWERS. A conversion can change the network quietly: an
operator mapped slightly differently, a layer's layout transposed. So this
also writes <stem>.check.json -- a few fixed test inputs and what TensorFlow
predicts for them. owg_live.py --check feeds the same inputs to the exported
file on the station and must get the same numbers back (within 1 mm) before
any of its wave heights mean anything.

Usage (Colab):
    !pip install -q tf2onnx onnx
    !python export_owg_onnx.py --model "/content/drive/MyDrive/Optical Waves/owg_c2_H_current_C"
"""

import os
import sys
import json
import argparse
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--model", required=True,
                    help="Model stem: reads <stem>.report.json and <stem>.weights.h5")
    ap.add_argument("--opset", type=int, default=13)
    ap.add_argument("--optimise", action="store_true",
                    help="Run tf2onnx's ONNX optimiser (needs ~10 GB of memory, more than "
                         "free Colab has; it only makes inference ~30%% faster)")
    args = ap.parse_args()
    # CPU only: the conversion gains nothing from the GPU, and TensorFlow's GPU
    # and XLA set-up costs memory Colab cannot spare.
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

    stem = args.model
    rep = json.loads(Path(stem + ".report.json").read_text())
    width = int(rep["img_size"])
    height = int(rep.get("img_height") or width)
    print(f"model             : {Path(stem).name}, {rep['arch']}, input {width}x{height}, "
          f"crop {rep.get('crop')}")

    import tensorflow as tf
    import tf2onnx
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from train_marconi_owg import build_model

    model = build_model(height, width, rep["arch"], 0.5)
    model.load_weights(stem + ".weights.h5")

    # Fixed test inputs, normalised like real images (mean 0, sd 1 each).
    rng = np.random.default_rng(0)
    X = rng.standard_normal((4, height, width, 1)).astype(np.float32)
    X = (X - X.mean(axis=(1, 2, 3), keepdims=True)) / X.std(axis=(1, 2, 3), keepdims=True)
    y_tf = model.predict(X, verbose=0).ravel()

    # Converted from a plain tf.function rather than from_keras, which
    # trips over Keras 3 models (Colab's TensorFlow 2.20).
    spec = (tf.TensorSpec((1, height, width, 1), tf.float32, name="image"),)

    @tf.function(input_signature=spec)
    def predict(image):
        return model(image, training=False)

    # tf2onnx's ONNX optimiser copies the whole graph, weights included, for
    # each of its passes: ~10 GB for Inception-ResNetV2, which crashed free
    # Colab (12.7 GB). Without it the file is equivalent (checked below and
    # by owg_live.py --check), just less fused.
    if not args.optimise:
        tf2onnx.convert.optimizer.optimize_graph = lambda graph, *a, **k: graph

    onnx_path = stem + ".onnx"
    tf2onnx.convert.from_function(predict, input_signature=spec, opset=args.opset,
                                  output_path=onnx_path)
    print(f"wrote             : {onnx_path} ({Path(onnx_path).stat().st_size / 1e6:.0f} MB)")

    try:
        import onnxruntime as ort
        sess = ort.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
        name = sess.get_inputs()[0].name
        y_onnx = np.array([sess.run(None, {name: X[i:i + 1]})[0].ravel()[0]
                           for i in range(len(X))])
        diff = float(np.abs(y_onnx - y_tf).max())
        print(f"ONNX vs TensorFlow: largest difference {diff * 1000:.2f} mm "
              f"({'OK' if diff < 1e-3 else 'TOO LARGE -- do not use this file'})")
    except ImportError:
        print("onnxruntime not installed here; the station's --check will compare instead")

    np.save(stem + ".check_inputs.npy", X)
    Path(stem + ".check.json").write_text(json.dumps(
        {"inputs": Path(stem).name + ".check_inputs.npy",
         "expected": [float(v) for v in y_tf]}, indent=1) + "\n")
    print(f"wrote             : {stem}.check.json and {Path(stem).name}.check_inputs.npy")
    print("Copy these four files to the station: .onnx, .report.json, .check.json, "
          ".check_inputs.npy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
