# This file is for validation of the GNN model on the test dataset.

import torch
from torch.cuda.amp import autocast
from torch_geometric.loader import DataLoader

from pathlib import Path

import time
from datetime import datetime
def val_loop_optimized(val_data, best_model=None, use_amp=True, device='cuda', model=None, criterion=None, batch_size=32,
                       exp = None, train_id = None,logs_dir = '/home/mriveraceron/glv-research/GNN-Logs'):

    # Section: Save lines
    val_lines = []
    start_time = time.time()
    now = datetime.now()
    formatted_time = now.strftime("%Y-%m-%d %H:%M:%S.%f")
    val_lines.append(
        f"{'-' * 40}\n"
        f"Starting validation\n"
        f"Timestamp: {formatted_time}\n"
        f"Experiment ID: {exp}\n"
        f"Version: {train_id}\n"
        f"The batch size for validation is: {batch_size}\n"
        f"{'-' * 40}\n"
    )

    # Create DataLoader for validation data
    val_loader = DataLoader(val_data,
        batch_size=batch_size,  # Can be larger than training since no gradients
        shuffle=False,  
        num_workers=4,
        pin_memory=True,  # Faster GPU transfer
        drop_last=False  # Keep last batch even if smaller than batch_size
    )
    # Load model weights if path is provided
    if best_model is not None:
        best_params = torch.load(best_model, weights_only=True)
        model.load_state_dict(best_params['model_state_dict'])
        model.to(device)

    model.eval()
    total_val_loss = 0
    num_samples = 0
    is_cuda = next(model.parameters()).is_cuda

    # Track memory usage
    if is_cuda:
        torch.cuda.reset_peak_memory_stats()
        initial_memory = torch.cuda.memory_allocated() / 1024**2  # MB

    try:
        with torch.no_grad():
            for data in val_loader:
                data = data.to(device)
                batch_size_actual = data.y.size(0)
                if use_amp and is_cuda:
                    with torch.autocast(device_type=device):
                        out = model(data)
                        loss = criterion(out, data.y)
                else:
                    out = model(data)
                    loss = criterion(out, data.y)
                
                # FIXED: Weight loss by actual batch size for correct averaging
                total_val_loss += loss.item() * batch_size_actual
                num_samples += batch_size_actual    
    except RuntimeError as e:
        if "out of memory" in str(e):
            error_msg = f"CUDA OOM error during validation. Try reducing batch_size from {batch_size}\n"
            val_lines.append(error_msg)
            print(error_msg)
            if is_cuda:
                torch.cuda.empty_cache()
            raise
        else:
            raise

    avg_val_loss = total_val_loss / num_samples
    elapsed_time = time.time() - start_time

    # Compute performance metrics
    metrics = {
        'loss': avg_val_loss,
        'num_samples': num_samples,
        'num_batches': len(val_loader),
        'time_seconds': elapsed_time,
        'samples_per_second': num_samples / elapsed_time
    }

    if is_cuda:
        peak_memory = torch.cuda.max_memory_allocated() / 1024**2  # MB
        metrics['peak_memory_mb'] = peak_memory
        metrics['memory_delta_mb'] = peak_memory - initial_memory

    # Section: Logging
    line = f'>> Validation Loss: {avg_val_loss:.6f}'
    print(line)
    val_lines.append(line + '\n')

    line = f'>> Samples Evaluated: {num_samples} (in {len(val_loader)} batches)'
    print(line)
    val_lines.append(line + '\n')

    line = f'>> Validation Time: {elapsed_time:.2f}s ({metrics["samples_per_second"]:.2f} samples/s)'
    print(line)
    val_lines.append(line + '\n')

    if is_cuda:
        mem_line = f'>> Peak Memory: {metrics["peak_memory_mb"]:.2f} MB (Δ{metrics["memory_delta_mb"]:.2f} MB)'
        print(mem_line)
        val_lines.append(mem_line + '\n')

    if is_cuda:
        line = f'>> Peak Memory: {metrics["peak_memory_mb"]:.2f} MB'
        print(line)
        val_lines.append(line + '\n')

    # Section: Save-logs
    if exp and train_id:
        try:
            log_path = Path(logs_dir) / exp / f'{train_id}_Validation.log'
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, 'w') as log_file:
                log_file.writelines(val_lines)
            print(f"Validation log saved to: {log_path}")
        except Exception as e:
            print(f"Warning: Could not save log file: {str(e)}")

    return metrics