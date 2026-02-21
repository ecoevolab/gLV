# --------------------------
# 24- January-2026
# This script defined the training function for GNN models.
# --------------------------
import torch
import time

def training_fn(model, device, batched_paths, loss_fn, optimizer, epochs=100):
    model.train()
    loss_history  =  []         # Loss at epoch
    total_elapsed = 0           # Running time
    best_epoch = 0              # Best epoch
    best_loss = float('inf')    # Best loss
    for iter in range(1, epochs+1):
        start = time.time()
        epoch_loss = 0
        for path in batched_paths:
            data_list = torch.load(path, weights_only=False)          
            for data in data_list:
                data = data.to(device)
                optimizer.zero_grad()
                out = model(data)
                loss = loss_fn(out, data.y)
                loss.backward()
                #torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item()   # Accumulate loss
        # Section: Best loss
        best_loss = float('inf') if iter == 1 else best_loss
        best_loss = min(best_loss, epoch_loss)
        best_epoch = iter if best_loss == epoch_loss else best_epoch
        # Append epoch loss to history
        loss_history.append(epoch_loss)
        elapsed = time.time() - start
        total_elapsed += elapsed
        # Print every 25 epochs
        if iter % 25 == 0:
            print(f"Epoch {iter}: Loss = {loss},  Elapsed time: {elapsed:.2f}")
    # Summary
    print(f'>> the total elapsed time with {epochs} epochs is {total_elapsed:.2f} seconds ( {total_elapsed/60:.2f} minutes)')      
    return  loss_history, best_loss, best_epoch, total_elapsed