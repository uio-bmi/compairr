/*
    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

class ThreadRunner
{
private:

  int64_t thread_count;

  pthread_attr_t attr;

  struct thread_s
  {
    int64_t t;
    void (*fun)(int64_t t);
    pthread_t pthread;
    pthread_mutex_t workmutex;
    pthread_cond_t workcond;
    int64_t work; /* 1: work available, 0: wait, -1: quit */
  } * thread_array;

  static void * worker(void * vp)
  {
    struct thread_s * tip = static_cast<struct thread_s *>(vp);

    pthread_mutex_lock(&tip->workmutex);

    /* loop until signalled to quit */
    while (tip->work >= 0)
      {
        /* wait for work available */
        if (tip->work == 0)
          pthread_cond_wait(&tip->workcond, &tip->workmutex);

        if (tip->work > 0)
          {
            (*tip->fun)(tip->t);
            tip->work = 0;
            pthread_cond_signal(&tip->workcond);
          }
      }

    pthread_mutex_unlock(&tip->workmutex);
    return nullptr;
  }

public:

  ThreadRunner(int t, void (*f)(int64_t t))
  {
    thread_count = t;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    /* allocate memory for thread data */
    thread_array = static_cast<struct thread_s *>
      (xmalloc(static_cast<uint64_t>(thread_count) * sizeof(struct thread_s)));

    /* init and create worker threads */
    for(int64_t i=0; i<thread_count; i++)
      {
        struct thread_s * tip = thread_array + i;
        tip->t = i;
        tip->work = 0;
        tip->fun = f;
        pthread_mutex_init(&tip->workmutex, nullptr);
        pthread_cond_init(&tip->workcond, nullptr);
        if (pthread_create(&tip->pthread,
                           &attr,
                           worker,
                           static_cast<void*>(thread_array + i)))
          fatal("Cannot create thread");
      }
  }

  ~ThreadRunner()
  {
    /* ask threads to quit */
    /* sleep until they have quit */
    /* destroy threads */
    /* finish and clean up worker threads */

    for(int64_t i=0; i<thread_count; i++)
      {
        struct thread_s * tip = thread_array + i;

        /* tell worker to quit */
        pthread_mutex_lock(&tip->workmutex);
        tip->work = -1;
        pthread_cond_signal(&tip->workcond);
        pthread_mutex_unlock(&tip->workmutex);

        /* wait for worker to quit */
        if (pthread_join(tip->pthread, nullptr))
          fatal("Cannot join thread");

        pthread_cond_destroy(&tip->workcond);
        pthread_mutex_destroy(&tip->workmutex);
      }

    xfree(thread_array);
    pthread_attr_destroy(&attr);
  }

  void run()
  {
    /* wake up threads */
    for(int64_t i=0; i<thread_count; i++)
      {
        struct thread_s * tip = thread_array + i;
        pthread_mutex_lock(&tip->workmutex);
        tip->work = 1;
        pthread_cond_signal(&tip->workcond);
        pthread_mutex_unlock(&tip->workmutex);
      }

    /* wait for threads to finish their work */
    for(int64_t i=0; i<thread_count; i++)
      {
        struct thread_s * tip = thread_array + i;
        pthread_mutex_lock(&tip->workmutex);
        while (tip->work > 0)
          pthread_cond_wait(&tip->workcond, &tip->workmutex);
        pthread_mutex_unlock(&tip->workmutex);
      }
  }
};

